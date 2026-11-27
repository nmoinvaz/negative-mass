# Evolving the phantom black universes along their growing mode (#217)

Code [evolve.py](evolve.py) (kernel with the potential term, data, extraction of the mode, departure of Y, finders, bisection, series, fits, checks; the potential and the member are in
[evolve_potential.py](evolve_potential.py)). Derivation [evolve_derive.py](evolve_derive.py), design [evolution_design.md](evolution_design.md) with the coordinator's rulings. Output of record
[evolve_output.txt](evolve_output.txt), assembled from negative-mass-raw `I217-mac/post`; the runs and the criteria (written before them) are in `I217-mac/runs` and `I217-mac/README.md`.
Status words as in AGENTS.md: every slope, A\*ₕ and finder result is *numerical evidence*. No sentence here says a result is new, and none names the threshold by a type.

## Answer

Along the family member + A × (growing mode), both c = −0.3 and c = −1.5 have a threshold A\*ₕ that converges at second order in h, and the lifetime slope equals 1/W of #212 in t (the static
observer's Killing time) to 1.6e-4 at c = −0.3 and 1.0e-4 at c = −1.5 at h = 0.0025, within fit errors of 0.04 % and 0.14 %. At c = −0.3 the two sides form two named outcomes (a trapped region and an
anti-trapped region); at c = −1.5 the collapse side is named by the finders and the other side by their absence to t = 50, with a late far anti-trapped cell. **The control c = −1 does not behave as
#212 says:** it shows a lifetime law with slope 1/κ = 1.5708 to 2.4e-4 and A\*ₕ at second order, so the control is NOT MET as found; the independent 1D master equation does decay at c = −1.
The one-sided pulse gives the same slope at c = −0.3 (to 3.6e-3) and does **not** at c = −1.5 (43 % below).

## What was measured

- **Member and family.** The static black universe of [BF06] with its potential V(φ), b = 1, phantom, in the tortoise null gauge (e^{−f} = A, r = R, φ = √2 arctan ρ), on the two rays through a corner at ρ_c (the throat at c = −0.3,
  ρ_c = 1 at c = −1.5 and −1). Data = member + A × δ̂, δ̂ the mode's shape in the code's own gauge from the untuned member's truncation-seeded departure at |Y| = 2e-4, A = Y at the corner node.
- **Time.** On the member t = (x⁺ + x⁻)/√2 is the Killing time (e^{−f} = A); the diagonal is the static observer at ρ_c; T is t of the diagonal cell, so the slope is compared with #212's W in t.
- **Lifetime and fate.** The first row whose diagonal cell has |Y| ≥ 2e-5, with Y = δr − (ρR/√2)δφ; the fate is the sign of Y there. This is a linear-regime measure: it tests the growth rate in the nonlinear code. It does not
  by itself show that generic data approach the member. (3e-3, the first ruling, is nonlinear at c = −1.5: the mode is mostly gauge, δf ≈ 10 to 20 × Y; with it A\*ₕ was −6.6e-06, −3.6e-06, −4.3e-06 and the series saturated near 3e-6, `dev_43` to `dev_45`.)
- **Fits.** As #40: bisection of the fate to 1e-14, A = A\* ∓ 10^−e (e = 5 to 11.5, step 0.25), T = c_side + s ln(1/|A − A\*|), common slope; primary window |A − A\*| ∈ [1e-7, 1e-5] (9 points per side), four more windows, joint and
  each side: 15 fits. The series below 1e-7 runs and is **not fitted**. Three error bars: the OLS SE (statistical); the fit spread, the envelope of s ± 2 SE over the 15 fits (window and side choice); and the resolution difference of
  the two finest h. "Inside the spread" means 1/W lies in the envelope.
- **Resolutions.** h = 0.02 (coarse member, reported), 0.01, 0.005, 0.0025; X = 40 (c = −0.3), 50 (c = −1.5, −1).

## Code checks (`evolve.py check`, 10 of 10; mutations fire exactly their declared checks)

| id | what | result |
|---|---|---|
| E0 | kernel with the potential flag off and on at c = 0 equals `dn.evolve`: static, collapse with excision, expansion, numba and numpy | max abs 0, dead masks equal |
| E1 | right-hand side against the equations with V (sympy), c = −0.3, −1.5, −1, 0 | 7.2e-16 (bar 1e-13) |
| E2 | numba against numpy kernel with V | 0 (bar 1e-12) |
| E3 | c = 0 rays with a pulse against `dn.ellis_rays` | r 1.7e-9, φ 5.5e-13, f 0 |
| E4 | static members c = 0, −0.3, −1.5, −1, t ≤ 8, h = 0.04, 0.02, 0.01 | orders 1.97 to 2.11 (c = −1.5: 1.99, 1.99) |
| E5 | growth rate of Y from the truncation seed against #212's W, |Y| ∈ [1e-4, 1e-3] | c = −0.3: −9.8e-3 (h = 0.01), −3.4e-3 (0.005); c = −1.5: −4.9e-3, −7.9e-4 (bar 1.5 %) |
| E6 | 1D master equation, no double-null code | W of #212 to ≤ 2e-6 at c = −2, −1.5, −1.1, −1.01, −1.001; c = −1: −0.037 (decays) |
| E7 | shape of the code's Y(l) against #212's Y(ρ(l)), l ∈ [−6, 6] | 4.8e-3 (c = −1.5), 2.7e-3 (−1.01), 2.6e-3 (−1.1), 5.3e-3 (−0.3, through the throat's arc); bar 1e-2 |
| PL, PX | planted series through the same `bisect`, `series`, `fit_set` (40 seeds, two W) | 1/W inside the spread 40 of 40, A\* recovered 40 of 40; 1/(√2 W) rejected 40 of 40 |

The three independent checks of the mode (ruling 1): W (E5, and E9 below), the shape against #212's Y(ρ(l)) (E7), and the 1D master equation (E6). E6 shares only #212's V̂; E7 shares #212's background and rhs.
Mutations (`th_coef`, `vterm_sign` (the wrong-sign-V mutation), `vterm_coef`, `Vp_sign`, `V_sign_f`, `numpy_pert`, `time_unit`, `pstar_off`): the fired sets equal the declared ones in the runs of record (`evolve_output.txt`).
In `analyze`, `time_unit` fires E9 at both members; `pstar_off` raised in E8 (A\* changed sign at c = −0.3, h = 0.0025), which counts as fired (README, A2).

## A\*ₕ and its order (E8: met at c = −0.3 and −1.5)

| member | h = 0.02 | 0.01 | 0.005 | 0.0025 | orders | Richardson (order 2) |
|---|---|---|---|---|---|---|
| c = −0.3 | −7.89399e-07 | −1.92364e-07 | −4.75546e-08 | −1.18496e-08 | 2.037, 2.016, 2.005 | +5.2e-11 |
| c = −1.5 | −2.85416e-06 | −7.04113e-07 | −1.77257e-07 | −4.52288e-08 | 2.019, 1.990, 1.971 | −1.2e-09 |
| c = −1 | −3.39873e-06 | −8.29464e-07 | −2.05469e-07 | −5.12506e-08 | 2.035, 2.013, 2.003 | +1.6e-10 |

The continuum threshold of a family through the static member is the member (A\* = 0): the Richardson values are zero to the extrapolation error.

## The lifetime slope against 1/W (E9: met at c = −0.3 and −1.5)

Primary joint slope ± OLS SE, relative deviation, fit spread, window [1e-7, 1e-5] (n = 18). 1/W = 1.03084 (c = −0.3), 2.02359 (c = −1.5); 1/κ = 1.57080 (c = −1, the control).

| member | h | slope ± OLS SE | relative | fit spread | collapse side / other side alone |
|---|---|---|---|---|---|
| c = −0.3 | 0.02 | 1.02119 ± 0.02766 | −9.4e-3 | [0.89018, 1.15279] | 1.12764 / 0.91473 |
| | 0.01 | 1.03040 ± 0.00699 | −4.3e-4 | [0.99928, 1.06438] | 1.05722 / 1.00358 |
| | 0.005 | 1.03122 ± 0.00178 | +3.6e-4 | [1.02099, 1.04029] | 1.03797 / 1.02446 |
| | 0.0025 | 1.03101 ± 0.00043 | +1.6e-4 | [1.02764, 1.03369] | 1.03255 / 1.02948 |
| c = −1.5 | 0.02 | 2.07676 ± 0.03950 | +2.6e-2 | [1.68234, 2.69352] | 2.04810 / 2.10543 |
| | 0.01 | 2.03377 ± 0.01579 | +5.0e-3 | [1.85799, 2.25398] | 2.00102 / 2.06653 |
| | 0.005 | 2.02528 ± 0.00706 | +8.3e-4 | [1.94590, 2.11056] | 2.00654 / 2.04401 |
| | 0.0025 | 2.02338 ± 0.00291 | −1.0e-4 | [1.99256, 2.05635] | 2.01473 / 2.03203 |
| c = −1 | 0.02 | 1.58460 ± 0.04160 | +8.8e-3 | [1.38153, 1.81361] | 1.74101 / 1.42819 |
| | 0.01 | 1.57375 ± 0.01020 | +1.9e-3 | [1.52318, 1.62689] | 1.61203 / 1.53546 |
| | 0.005 | 1.57170 ± 0.00218 | +5.8e-4 | [1.55983, 1.58392] | 1.57948 / 1.56392 |
| | 0.0025 | 1.57042 ± 0.00045 | −2.4e-4 | [1.56589, 1.57250] | 1.57190 / 1.56894 |

E9 at the two finest h (0.005, 0.0025): each primary slope within 2 OLS SE of 1/W, 1/W inside each fit spread, and the two slopes consistent (agreement +2.0e-4 at c = −0.3, +9.4e-4 at c = −1.5). The error shrinks with h
at every member. The three error bars are different things: the OLS SE assumes independent residuals, the fit spread covers window and side choice, the resolution difference compares the two finest h.

**Floors** (reported, the series below the window is not fitted): the largest |A − A\*| below the window at which T leaves the window's log law by more than 0.25.

| member | h = 0.02 | 0.01 | 0.005 | 0.0025 | T at e = 11.5, both sides, h = 0.0025 |
|---|---|---|---|---|---|
| c = −0.3 | 5.6e-8 | 1e-11 | none | none | 16.13, 16.17 |
| c = −1.5 | 5.6e-8 | 5.6e-8 | 3.2e-8 | 1.8e-8 | 18.72, 20.44 |
| c = −1 | 1.8e-8 | 1.8e-8 | 5.6e-9 | 1.8e-9 | 17.30, 18.82 |

At c = −1.5 the floor sits at |A − A\*| ≈ 3e-8 to 1.8e-8 and T at e = 11.5 rises by 1.7 to 2.0 per halving of h (13.3/14.8, 15.0/16.6, 16.8/18.5, 18.7/20.4), which is a seed that scales as about h^1.3 (ln 4/W = 2.8 would be h²);
its cause is not found. At c = −0.3 the series is clean to 1e-11.5 for h ≤ 0.005.

## The outcomes (naming runs, h = 0.005; finders raw = `dn.cell_class`, weighted = signs and 2e^f r₊r₋ ≥ 0.02)

| member | side | first trapped cell | first anti-trapped cell | the diagonal later |
|---|---|---|---|---|
| c = −0.3 | A < A\* (collapse) | weighted l = +0.16, r = 0.88, t = 7.1, 9.5, 11.8 for |A − A\*| = 1e-4, 1e-5, 1e-6 (raw: l = +0.02, r = 0.99, t = 5.1, 7.5, 9.8) | weighted l = −1.34, r = 1.10, t = 9.3, 11.6, 14.0 | excised |
| | A > A\* (expansion) | none (raw and weighted) | weighted l = −0.14, r = 1.13, t = 7.6, 10.0, 12.4 (raw: l = −0.02, r = 1.015, t = 5.2, 7.6, 10.0) | r = 1.6, 5.9, 10.7, 15.6, 20.6 at t = 10, 20, 30, 40, 50; r₊ → 0.70, r₋ → −7e-10 at t = 40 |
| c = −1.5 | A < A\* | weighted l = +0.90 to +1.45, r ≈ 1.0, t = 13.2, 17.9, 21.9 (raw: l = −0.02 to +0.28, r = 1.04, t = 12.3, 16.7, 20.7) | weighted l = +0.57 to +1.0, r = 0.73 to 0.81, t = 13.3, 17.6, 21.7 | excised |
| | A > A\* | none (raw and weighted) | raw none; weighted only at t = 60.2, 62.0, 64.5 (l = +6.5, +4.1, +2.0, r = 25.7, 23.3, 21.4) | r = 3.8, 8.1, 12.8, 17.6 at t = 20, 30, 40, 50 (|A − A\*| = 1e-4); r₋ → −5e-10 |
| c = −1 | A < A\* | weighted l = −0.5 to −0.3, r = 0.92, t = 10.6, 14.1, 17.6 | weighted l = −0.6 to −0.34, r = 0.69 to 0.73, t = 11.3, 14.8, 18.2 | excised |
| | A > A\* | none | weighted only at t = 50.9, 59.2, 61.3 (l = +3.7, +8.7, +6.8, r = 21.2, 26.1, 24.4) | as c = −1.5 |

- **Both outcomes named by the finders: met at c = −0.3**; at c = −1.5 (and c = −1) the collapse side is named by the finders and the other side by their absence through t ≈ 50 with a late far anti-trapped cell, so the second outcome is
  named by what the finders do not see.
- The collapse side forms a trapped and an anti-trapped region in the same run at every member. The expansion side at c = −0.3 forms an anti-trapped region and no trapped one.
- On the A > A\* side at c = −1.5 the diagonal's Misner–Sharp mass is 0.45, 0.38, 0.36 at t = 20, 30, 40 (r = 3.8, 8.1, 12.8; |A − A\*| = 1e-4) against the member's ADM mass m = 0.318, with 2M/r = 0.057 at t = 40. At c = −0.3
  (m = 0.064) the expansion side has M = 1.9, 1.9 at t = 20, 30 (r = 5.9, 10.7), 30 times m. M and r₋ at t = 50 are roundoff (e^f ≈ 1e9 times r₋ ≈ 1e-11); I quote nothing from t = 50.
  These are values, not an identification of the final state.
- Raw against weighted: at |A − A\*| = 1e-8 (h = 0.02, development, `dev_50`) the raw finder is silent on the expansion side of c = −0.3 and the weighted one fires at t = 17.3; in the naming runs the raw finders are silent for the expansion side
  at c = −1.5 and c = −1 at every |A − A\*|.

## The control: c = −1 (E11: NOT MET as found)

Declared before the runs (README): #212 has no growing mode at c = −1, so along the family through this member there is no divergent lifetime at A = 0. E11 is MET if no lifetime law is found, where "found" is: the primary slope's OLS SE ≤ 10 % at the
two finest h, the two slopes consistent, A\*ₕ at second order. **Found:** the slope is 1.57170 ± 0.00218 and 1.57042 ± 0.00045 (h = 0.005, 0.0025), 1/κ = 1.57080 to 5.8e-4 and 2.4e-4, A\*ₕ orders 2.035, 2.013, 2.003.
The constraints stay at truncation level while Y grows (`dev_76`), the growth is the same in a horizon-regular gauge (`dev_74`), and the independent 1D master equation decays (E6, −0.037; to t = 60 in `dev_30`). This is reported, not fixed.
**The open question** (follow-up #232): at c = −1 the double-null problem and #212's master equation disagree. What D1 and D1b add, and what they do not, is next.

## The one-sided pulse (family O, XC)

| member | h | slope ± OLS SE | mode family at the same h | relative difference |
|---|---|---|---|---|
| c = −0.3 (+ ray) | 0.01 | 1.04681 ± 0.00226 | 1.03040 | +1.6e-2 |
| | 0.005 | 1.03494 ± 0.00095 | 1.03122 | +3.6e-3 (XC met, bar 5 %) |
| c = −1.5 (− ray) | 0.01 | 1.14137 ± 0.03026 | 2.03377 | −4.4e-1 |
| | 0.005 | 1.15004 ± 0.03718 | 2.02528 | −4.3e-1 (XC **not met**) |

At c = −0.3 the pulse's collapse time follows the same law (T = 13.4 → 28.6 over e = 5 → 11.5, a steady 1.2 per half decade). At c = −1.5 it does not: T(e) = 24.0, 25.9, 27.5, 28.8, 30.0, 31.1, 32.2, 33.3, 34.5, 36.3 for e = 5 to 9.5
(increments 1.9 falling to 1.1 per half decade, then 48.5 at e = 10): the first trapped cell comes long after the departure and its delay depends on A − A\*, so a finder-event slope is not the departure slope there.
This is a measurement of a different quantity at c = −1.5, not a measurement of the slope's dependence on the family; the cross-check is not met and is reported as found.

## The response to a compact pulse (reported diagnostic, ruling 2)

D1 (declared): the odd part of Y in ±A_p = 1e-6 on the + ray (centre 3, width 0.5), h = 0.01, X = 40, to the end of the grid (background: the untuned member).

| member | rate [10,16] / [16,22] / [22,28] | Y/A_p late (median of the last 8 time units) | saturated from t | ratio to the response at 10 A_p, t = 54 | δM_diag/A_p, δM_edge/A_p at t = 34 (edge l) |
|---|---|---|---|---|---|
| c = −1.5 | +0.495 / +0.415 / −0.096 | 8.46e3 | 41.6 | 0.998 | −1.5e3, −7.0e2 (+11.3) |
| c = −1.01 | +0.405 / +1.089 / +0.136 | 9.93e3 | 34.1 | 1.004 | −1.0e3, −6.8e2 (+11.3) |
| c = −1 | +0.787 / +0.902 / +0.130 | 9.97e3 | 33.9 | 1.004 | −1.0e3, −6.8e2 (+11.3) |
| c = −0.99 | +0.678 / +0.812 / +0.125 | 1.00e4 | 33.8 | 1.004 | −1.0e3, −6.8e2 (+11.3) |
| c = −0.3 | +0.776 / +0.081 / +0.010 | 5.77e4 | 21.6 | (the 10 A_p run excises cells) | +2.8e3, +3.3e3 (+11.3) |

**Plainly: c = −1.5 and c = −1.01, which have #212 modes, saturate in this diagnostic, and the response stays linear in A_p (to 0.4 %) while it does.** That is not a statement about the member's linear problem: the untuned background has left the
member before the saturation. The time at which the untuned |Y| on the diagonal reaches 2e-5, 1e-3, 1e-2, 0.1 (`bgdep`): c = −1.5: 6.8, 14.8, 19.4, 23.4; c = −1.01: 5.0, 11.3, 14.9, 18.7; c = −1: 5.0, 11.2, 14.8, 18.6; c = −0.99: 5.0, 11.1, 14.7, 18.5;
c = −0.3: 4.7, 8.9, 11.2, 13.8. The odd response is linearized about that departing solution. The mass at the right edge changes by about −6.8e2 A_p at t = 34 at the four members with ρ_c = 1, and at the diagonal by −1.0e3 to −1.5e3 A_p, the same sign and no decay,
which is how a shift of the state would look and not how a bounded mode would (δM_edge → 0), but it is the departed solution that shifts.

**D1b** (added after D1, README A1): the background is the mode family at A = A\*ₕ (h = 0.01), so it stays at the member until its floor. The odd response then grows at (t = 6 to 10, 10 to 14, 14 to 18, 18 to 22):
c = −0.3: 0.968, 0.970, 1.055 (W = 0.97008); c = −1.5: 0.537, 0.501, 0.517, 0.725 (W = 0.49417), and it leaves at the floor (t ≈ 18 to 22); no saturation before that. At c = −1.01, −1, −0.99 the rates over the same windows are
−0.20, 0.86, 1.10, 2.17; 0.08, 0.71, 1.35, 2.23; 0.31, 0.42, 1.72, 2.29: not clean. Their floor departure (T ≈ 14 at h = 0.01) overlaps the windows, and the response may also still carry the pulse's slow tail at the start; I have not separated the two.
So: **at c = −0.3 and c = −1.5 the response about the member follows #212's W up to the floor and does not saturate in that range; for the three members near c = −1 this diagnostic is inconclusive**, and the late-time saturation of D1 is the background's departure.
The statement that the double-null linear problem disagrees with the master equation "at late times at every member" is therefore not supported; the c = −1 disagreement rests on the untuned growth at κ and the lifetime law (above), not on D1.

## Acceptance of #217

| criterion | verdict | numbers |
|---|---|---|
| both outcomes found and named by the horizon and anti-trapped-surface finders | met at c = −0.3; at c = −1.5 the collapse side by the finders, the other by their absence and a late far cell | naming table |
| slope equals 1/W within a stated fit error, at two resolutions | met at c = −0.3 and c = −1.5 | +3.6e-4, +1.6e-4 (c = −0.3); +8.3e-4, −1.0e-4 (c = −1.5), h = 0.005, 0.0025; fit errors 0.17 % / 0.04 % and 0.35 % / 0.14 % |
| control c = −1 shows no such threshold | **NOT MET** as found | lifetime law with slope 1/κ to 2.4e-4, A\*ₕ second order |
| potential term passes a planted static check | met | E0 to E4, mutations fire |
| no sentence calls the threshold by a type; only outcomes and slopes | met | |
| cross-check: the one-sided pulse gives the same slope (XC) | met at c = −0.3, not met at c = −1.5 | +3.6e-3; −43 % |

## What ran

16 jobs (4 resolutions at each of three members, mode family; two resolutions at two members, family O), exit codes 0; times 7 s to 1064 s (the h = 0.0025 jobs 449 s, 819 s, 700 s at c = −0.3, −1.5, −1). `evolve.py check` rc 0 (10 of 10);
`analyze` rc 1 because E11 and XC[c = −1.5] are not met, as declared; eight mutated checks rc 1; naming and D1 rc 0. Two processes, nice 10, one thread each. Development: `dev_01` to `dev_122` in the raw folder.

## What this does not show

- A lifetime slope that equals 1/W in a family through the static member tests the growth rate in a nonlinear code. It does not show that generic data approach the member (a family not through it would, as #40's G).
- The mode is extracted from the code; it is checked against #212 (E5, E6, E7) and not derived from #212's gauge. Only the radial mode of two members and a control, D = 4, one potential. Nothing about polar modes (#218) or the sourced black bounce (#223).
- The horizon is never reached (tortoise gauge); the interior is not evolved. What state the side without a finder event ends in is not identified; late masses are roundoff-limited beyond t ≈ 40.
- At double precision the member cannot be followed past the floor (T ≈ 16 to 20): beyond it the seed (roundoff and truncation) departs whatever A is. Statements about "late times" at the member are limited to that range.

## Recommended tasks

1. **#232, the double-null linear response against #212's master equation at c = −1.** Proposal for that issue. *What:* (a) cross the horizon: the Kruskal data of the member in both null coordinates including the interior, and run c = −1.01, −1, −0.99 across it
   (`evolve_dev_reg.py` has the two-step-size kernel and the exterior part); (b) a direct linear double-null solver: linearize the three evolution equations about the member in the tortoise gauge, evolve constraint-satisfying compact data, extract Y, and compare it with the 1D
   solution of #212's equation at c = −1.5, −1.01, −1, −0.99 over t ≤ 14, where the background of the nonlinear code is still the member. *Why:* the nonlinear code has a lifetime law at c = −1 that the master equation does not; (b) separates the linear problem from the code. *Candidate for the reader:*
   #212's variable u = Y/ρ and its δφ = 0 gauge at the member whose horizon sits at the throat, ρ_h = 0, where the horizon, the throat's wall 2A²/ρ² and R′ = 0 meet and δλ = (Y′ − A′Y/2A)/R′ divides by R′; the decaying solution of u may not be the regular perturbation there.
   *Acceptance (each can fail):* the linear solver reproduces the 1D growth rates at c = −1.5 and −1.01 to 1e-3 (a control that fires if the solver is wrong); at c = −1 it either grows at κ (then the master equation's no-mode count needs a separate look) or decays (then the nonlinear code's growth is a discretization or gauge effect, and the
   regular-gauge result of `dev_74` needs the interior). *Depends on:* #217, #212. The late-time clause of the first statement of the question is withdrawn (D1b).
2. **The c = −1.5 floor.** *What:* find why a seed ≈ h^1.3 that A cannot tune departs at T ≈ 15 to 20 (|A − A\*| ≈ 2e-8 to 3e-8). *Why:* it limits the window and the dynamic range. *Acceptance:* a run in extended precision (or with the data perturbed at 1e-16) moves the floor
   by the amount the seed predicts, and the h exponent is explained by a named term. *Depends on:* #217, #40's roundoff task.
3. **A family-independence check at c = −1.5.** *What:* a one-sided pulse placed so that its response is gone before the departure is measured (or the pulse's departure by the same Y₀), at c = −1.5 and −1. *Why:* XC could not test family independence there (finder events lag by an A-dependent delay).
   *Acceptance:* the pulse's departure slope at h = 0.005, 0.0025 equals the mode family's within 5 %; a mutation that scales T by √2 fails it.
4. **The mode by the analytic route.** *What:* build #212's Y(ρ) on the whole static region (arc at c = −0.3), the gauge transformation to ray data with the δf = 0 ray gauge, and use it as a third family. *Why:* an independent δ̂. *Acceptance:* its A\*ₕ is second order and its slope equals the extracted family's within the fit error.
5. **The state on the A > A\* side at c = −1.5 and c = −1.** *What:* identify it (mass against r at t ≤ 40, r₊ → 1/√2, r₋ → 0⁻) with a run that does not lose r₋ to roundoff (a rescaled gauge). *Acceptance:* its M(r) matches a named static solution of the same potential to 1e-3, or it is shown not to.

## Files and reproduction

`evolve.py` (this directory, `~/Source/negative-mass/.venv/bin/python`, from `notes/black_universe`), `evolve_potential.py`, `evolve_derive.py` (+ `evolve_derive_output.txt`), `evolve_output.txt`, `evolve_dev.py`, `evolve_dev_reg.py`, `evolve_dev_runs.py`, `evolve_dev_runs2.py`,
`evolve_master1d.py`, `evolution_design.md`, this note. Raw: negative-mass-raw `I217-mac` (README with the criteria, `runs/`, `post/`, `dev/`, `scripts/`).

```
python evolve.py check                       # E0 to E7 and the planted series, 20 s
python evolve.py check --mutate=vterm_sign   # a mutation; must fire E1, E4, E5, E7
python evolve.py job --c=-1.5 --h=0.005 --out=DIR              # about 3 minutes; --family=pulse for family O
python evolve.py analyze --raw=I217-mac/runs                   # E8, E9, XC, E11 and the tables
python evolve.py name --c=-0.3 --h=0.005 --raw=I217-mac/runs --out=DIR
python evolve.py resp --c=-1.5 [--tuned=I217-mac/runs/c-1.5_h0.01]
```

The runs of record took about 25 minutes per lane on two processes (the post-processing 8 minutes).

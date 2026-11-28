# The settled black-hole mass at the Ellis threshold, and the masses 0.30 a and 0.22 a (#41)

Code: [settled_mass.py](settled_mass.py) (imports #39's `dn.py` and #40's `threshold.py` unchanged). Output of record:
[settled_mass_output.txt](settled_mass_output.txt). Raw output, with the criteria written before the runs of record:
negative-mass-raw `I41-rpi5` (README, `runs/`, `mut/`, `dev/`). Status words as in AGENTS.md: every mass below is *numerical evidence*;
what is not shown is marked *conjecture*. No sentence here says a result is new. Units: a = 1 is the throat radius of the initial static
wormhole, 8πG = 1, D = 4, and M = r(1 − (∇r)²)/2, which on a marginally trapped sphere (r₊ = 0) is r/2.

## Answer

In the symmetric ghost-pulse family S at D = 4 the Misner–Sharp mass on the apparent horizon falls from 0.494 a at the first trapped cell and
settles at **M_AH = 0.22308 a (r_AH = 0.4462 a)**, with a spread of 3e-6 over resolution, window and δ, to the declared limit e^{−f} = 1.57e6 on
the horizon, 9.4 time units after formation; it is the same for p*ₕ − δ with δ = 1e-8 and 1e-6 (and for δ up to 1e-3). This is [GGS09b]'s
r_AH = 0.45 (m_AH ≈ 0.22). [SH02]'s 0.30 a is not reproduced, even with SH02's own initial data run through the same code; [ST17]'s number
is the formation value. **The discrepancy is resolved for [GGS09b] and [ST17] and remains for [SH02]**: the units and definitions agree,
and the readings differ in time.

## A. What each paper calls the threshold black-hole mass

All three were read from arXiv (PDF; TeX source for [GGS09b] and [ST17]); hashes and fetch log in the raw folder (`papers.sha256`,
`papers.fetch.log`); the files are kept outside the repository.

| | [SH02] (gr-qc/0205041) | [GGS09b] (0806.1370) | [ST17] (1706.02070) |
|---|---|---|---|
| quoted | "M ∼ 0.30a" (3.6); Table I: 0.296 at E₀/a = 4.39e-5 | "m_AH ≈ 0.22" (Sec. IV C); "r_AH = 0.51 for ε_c = 0.07 to r_AH = 0.45 for ε_c = 0" (Sec. VI) | E_f/a₀ = 3.14, 6.28, 9.87 for n = 4, 5, 6 (Table I, Sec. III B), "does not depend on the injected energy" |
| definition | M = E(∞, x_H⁻) (3.4); E = ½r[1 − g⁻¹(dr, dr)] = ½r + e^f r ∂₊r ∂₋r (2.36), Sec. II D | M = (r/2)(1 − g̃^{ab}∇ₐr∇_b r) (22), Sec. III C, on the apparent horizon, "an outermost marginally trapped surface" (Sec. III B) | E_n of (2.47), Sec. II G: (n−2)A_{n−2}/(2κ²Ω){k + (2e^f/(n−2)²)ϑ₊ϑ₋ + …}, A_{n−2} as printed |
| unit | G = 1 (below 2.33); a = throat radius; E(throat) = a/2 (2.39) | G = 1 (M = r/2 on r = 2M); b = 1 = throat radius (Sec. III) | a₀ = throat radius; A_{n−2}, κ as printed |
| read | at x⁻ = x_H⁻ for x⁺ = 12, 16, 20 and ∞ (Fig. 7) | r_AH(t) settles by t ≈ 10 to 14 (Fig. 7), "final" | at the end of the grid, x⁺ = 5 (Table I caption) |
| data | one-sided ghost pulse on Σ₊, ℘₊ = a/√(2a² + (x⁺)²) + c_a exp(−c_b(l − c_c)²), c_b = 3, 6, 9, c_c = 3 (3.1); E₀ = E(20, 0) − static (3.3) | time-symmetric, c̄ = c̄_static + ε_c exp(−(x − x_c)²/σ_c²) (19), σ_c = 0.5, x_c = 0 | normal pulse π₊ = c₁ exp(−c₂(x⁺ − c₃)²), c₂ = 16, c₃ = 0.7; GR rows c₁ = 0.25, 0.5 |
| in this unit | M = 0.296 to 0.30 a (r = 0.59 to 0.60 a) | r_AH = 0.45 gives M = 0.225 a; 0.51 gives 0.255 a | E_f/E_throat = 3.14/π = 1.00: M = 0.50 a, r = a₀ |

Conversions:

- **[SH02].** (2.36) is `dn.mass_nodes` for n = 2. Their G = 1 against this code's 8πG = 1 changes the field normalization (φ_code = √2 φ_SH, from
  R_ab = ∓2∂φ∂φ against ∓∂φ∂φ) and not E, which is geometric. No factor.
- **[GGS09b].** The same geometric M, G = 1. The text prints m_AH ≈ 0.22 without defining m_AH. At a marginally trapped sphere (22) gives M = r/2, and Sec. IV C
  evaluates K = 2M/r on the horizon and finds it tends to 1. So m_AH = r_AH/2 and 0.45/2 = 0.225, against 0.22 printed. (An inference from (22); the text
  never writes m_AH = r_AH/2.)
- **[ST17].** E_n is normalized to π a₀ (n = 4) at the static throat, so 3.14 is the static throat's own energy: the static part of E_i at x⁺ = 5 is
  E_i − ΔE = 0.85 = 0.2722 π (0.2722 = a/r at l = 5/√2), for both rows, and E_f = 3.14. E_f/E_throat = 1.00 at n = 4, 5, 6. The "final mass" is the mass of a
  horizon at r = a₀, read 2 to 3 a₀ after x_H⁻ (2.1 to 2.9): the formation value, and that is why it "does not depend on the injected energy". Discrepancy in the
  paper, not fixed: the printed A₄ = π²/Ω⁴ gives E_f/a₀ = 2π² = 19.7 at n = 6, the table has π² = 9.87, which needs A₄ = π²/(2Ω⁴) (the volume of the unit 4-ball).

Where the printed times differ: [ST17] reads 2 to 3 a₀ after x_H⁻, [SH02] 5 to 17 time units after formation, [GGS09b] "late". The next section shows M_AH is still falling 3 time units
after formation.

## B. The settled mass in this code

**Method** (not in #39 or #40). A uniform grid cannot follow the horizon: e^{−f} grows exponentially along it (1.0 at formation, 4.1 at +2.5 time units, 100 at +4, 1.6e6 at +9.4),
the horizon is a layer of width about r e^{f} in x⁻, and #39's cell update loses it (F_cut = 25). The grid is a chain of levels:

```
 x⁻ ▲   level 2 (width h/4, Nc cells, starts at v₂)        box k+1: lower ray = a row of box k
    │        ┌───────┐                                      left ray  = box k's column at v_{k+1},
    │     ┌──┴───────┴──┐                                              cubic-interpolated to the finer nodes
    │     │ level 1 (h/2, starts at v₁)                     v_{k+1}: first column where e^{−f} on the horizon
    │  ───┴───╲ apparent horizon (r₊ = 0)                            reaches F₁ 2^k, so e^{−f}h₋ is the same at the
    │ level 0: uniform (h, h), #39's grid, bit-identical (PL3)       end of every level
    └─────────────────────────────────────────────────────▶ x⁺
```

Characteristic evolution has no dependence on later x⁻ or x⁺, so each box is exact in its domain of dependence. F₁ = 1.5, 20 levels, the declared limit is
e^{−f} ≤ F_max = 1 572 864 on the horizon; level 1 is 2.4 wide in x⁻ (Nc = 4.8/h cells). The horizon is the + to − zero of r₊ in each cell column, linearly
interpolated, and M_AH = r/2 there; the series starts at #40's first trapped cell (r₊, r₋ < −1e-2). Parameters and how they were fixed in development:
raw README.

**Controls** (all in `settled_mass_output.txt`; 13 of 13 pass in the runs of record, and each of 15 mutations fires exactly its declared set, in
development and in the runs of record).

| id | what | result |
|---|---|---|
| PL1 | exact Schwarzschild, r_h = 0.6, in the stretching gauge, 20 levels | max \|M − 0.3\| = 1.1e-4, 2.8e-5, 6.8e-6, 1.7e-6 at h = 0.02 … 0.0025; orders 2.00, 2.02, 2.02 |
| PL2 | exact ghost solution, moving contour r₊ = 0, two forced zoom levels | max error of r 1.7e-4, 4.3e-5, 1.1e-5, 2.7e-6; orders 2.00, 2.00, 2.00 |
| PL3 | level 0 against `dn.evolve` | bit-identical nodes and dead flags, same first trapped cell (385, 385) |
| AN | chain against a zoom-free uniform anisotropic grid (h, h/32), e^{−f} < 40 | max \|ΔM\| = 2.0e-4 (bar 3e-4; each grid at its own p*ₕ) |
| RC1–RC3 | p*ₕ on this machine, fates, first trapped cell, formation value | p*ₕ within 5e-14 of #40's (dh610), cells (522, 522), (1045, 1045), (2090, 2090) equal #40's; M_AH at formation within 1.5e-4 of #40's M_hor |
| RC4 | first law ΔM = ∫(r²/2)e^f r₋φ₊² dx⁺ on the horizon, to e^{−f} = 20 | relative difference 3.9e-4 |
| RC5 | every run reaches the declared limit, ≥ Nc/8 cells from the box edges | F_end 1.53e6 to 1.57e6 |
| RC6, RC6b | second order of M_AH(x⁺), three windows, both δ | orders 1.89 to 2.18 at h = 0.02, 0.01, 0.005; 1.94 to 2.41 at 0.01, 0.005, 0.0025 |
| RC7, RC8 | independent of δ; settled at the limit | plateau means at δ = 1e-6, 1e-8 differ by 1.7e-5 (h = 0.02) to 7e-8 (h = 0.0025); change between the last two windows ≤ 1.7e-5; flux 2.6e-7 |

**M_AH(t)** (finest h = 0.0025, δ = 1e-8; x⁺_f = 10.4488, t_f = 14.763; δ = 1e-6 gives the same M at the same x⁺ − x⁺_f to 1e-5, with t_f = 10.867):

| x⁺ − x⁺_f | t | r_AH | M_AH | e^{−f} on the horizon |
|---|---|---|---|---|
| 0 | 14.763 | 0.98848 | 0.49424 | 1.02 |
| 1.0 | 16.121 | 0.94147 | 0.47073 | 1.13 |
| 2.0 | 17.286 | 0.75600 | 0.37800 | 1.94 |
| 2.5 | 17.728 | 0.60630 | 0.30315 | 4.07 |
| 3.0 | 18.103 | 0.49997 | 0.24998 | 12.0 |
| 3.5 | 18.459 | 0.45860 | 0.22930 | 36 |
| 4.0 | 18.813 | 0.44863 | 0.22432 | 100 |
| 4.5 | 19.167 | 0.44740 | 0.22370 | 260 |
| 5.0 | 19.520 | 0.44735 | 0.22367 | 653 |
| 6.0 | 20.228 | 0.44669 | 0.22335 | 3.95e3 |
| 7.0 | 20.935 | 0.44625 | 0.22313 | 2.31e4 |
| 8.0 | 21.642 | 0.44617 | 0.22308 | 1.31e5 |
| 9.4 | 22.632 | 0.44617 | 0.22308 | 1.44e6 |

M_AH is monotone decreasing (the first law: ghost flux, r₋ < 0). It crosses 0.30 at x⁺ − x⁺_f = 2.5 and is within 0.6 % of its settled value from 4.0. Between 4.7 and 8 the field φ₊ on the
horizon changes sign and a second, smaller infall of negative energy (flux up to 4e-4) lowers M by 6e-4; after 8.5 the flux is below 1e-6. With the limit at e^{−f} = 6144 (x⁺ − x⁺_f = 6) the flux had not died; that
is why the limit is 1.57e6 (raw README, "Changes made during development").

**Settled value and its error bars** (kept separate):

| kind | what | size |
|---|---|---|
| value | plateau mean over the last 1.0 of x⁺, h = 0.0025 | 0.223083 (δ = 1e-8), 0.223083 (1e-6); r_AH = 0.44617 |
| resolution | plateau means at h = 0.02, 0.01, 0.005, 0.0025 (δ = 1e-8) | 0.223155, 0.223102, 0.223085, 0.223083; Richardson (second order, three finest) 0.223082; finest pair differs by 2.1e-6 |
| window | range of M over the last 1.0, finest h; M at the limit against the plateau | 2.6e-6; 0.223085 against 0.223083 |
| δ | plateau at δ = 1e-6 against 1e-8, finest h | 7e-8 (1.7e-5 at h = 0.02) |
| limit | mean horizon flux over the last 0.5 | 2.6e-7 per unit x⁺; no extrapolation is made, and the earlier second infall shows the flux can return, so "settled" is "to e^{−f} = 1.57e6" |

M at the limit against the distance from threshold (family S, h = 0.01, `scan`): 0.2231 for δ = 1e-8, 1e-7, …, 1e-3 (A = p*ₕ − δ); 0.2228 at 1e-2, 0.2294 at 3e-2, 0.2332 at 1e-1. E₀ of
(3.3) is not the fate-fixing quantity for S (a Gaussian in φ has E₀ ∝ +A at first order, and the threshold is the symmetric combination), so the scan is not set against SH02's E₀.

## C. The comparison

| | printed | here, same unit | verdict |
|---|---|---|---|
| [GGS09b] | r_AH = 0.45 (ε_c → 0), m_AH ≈ 0.22 | r_AH = 0.4462, M = 0.22308 | agrees to their two digits (r 0.9 % below 0.45); same threshold (the static wormhole) |
| [ST17] | E_f/a₀ = 3.14 = E_throat, read at x⁺ = 5 | M = 0.494 at formation, 0.47 at +1 | the formation value, not a settled mass |
| [SH02] | M ∼ 0.30 a; 0.296, 0.310, 0.325, 0.420, 0.38, 0.35 (Table I) | 0.2232, 0.2237, 0.2279, 0.2603, 0.2507, 0.2462 with SH02's own data (h = 0.005) | **not reproduced**: 0.07 to 0.16 a lower |

**SH02's data through the same chain** (`sh02`, h = 0.01 and 0.005, which agree to 1e-4 in M; Gaussian in l = x⁺/√2 as printed in (3.1)): E₀ = 4.23e-5, 4.22e-4, 4.16e-3, 3.57e-2, 2.51e-2, 2.05e-2 against
their 4.39e-5, 4.38e-4, 4.32e-3, 3.70e-2, 2.60e-2, 2.12e-2 (all 3.4 to 3.8 % below; the pulse is read right). The horizon coordinate x⁻_H (h = 0.005) is 8.08, 6.92, 5.59, 4.27, 4.61, 4.78
against 8.42, 6.98, 5.64, 4.46, 4.82, 4.90. The masses rise with E₀ as theirs do, but sit lower. Reading the Gaussian at x⁺ = c_c instead of l = c_c gives E₀ 17 % low and was rejected on that.

**What could cause the difference**, with what was tested:

- *Units and definitions.* No. E is the same geometric quantity, G = 1 in both, a is the throat radius in both ([SH02] (2.36), (2.39); this code's `mass_nodes`).
- *The family.* No, to the extent tested. S gives 0.2231 over δ from 1e-8 to 1e-3, and SH02's own one-sided data give 0.2232 at c_a = −1e-4 (E₀ = 4.2e-5). Larger E₀ raises M (0.26 at E₀ = 3.6e-2)
  here and in their table, but their small-E₀ values (0.296 to 0.325) are above the plateau of this code at the same E₀.
- *The time of reading.* This is the one difference the three papers have, and M_AH(t) above shows its size. The e^{−f} at which this code's M_AH equals SH02's printed M is 4.5, 3.8, 3.4, 2.1, 2.7, 3.4
  (1.5 to 2.5 time units after formation) for the six rows: their numbers equal a horizon mass read about two time units after formation, the kind of reading [ST17] makes at x⁺ = 5. SH02 say they read at x⁺ = 12, 16,
  20, 5 to 17 time units after formation. **Conjecture:** at those x⁺ their uniform grid (9601 points over x⁺ = [0, 20], h = 2.1e-3) cannot resolve a horizon with e^{−f} of 1e3 or more (the cell
  changes r₊ by e^{−f}h/2r ≈ 1 or more), and the E they read there is not a resolved horizon's. Test: this code's kernel on one uniform grid of that size (`sh02_uniform`) loses the horizon at
  e^{−f} = 101 to 138 (x⁺ = 7.3 to 10.6, before their x⁺ = 12) and the last mass it finds there is 0.229, 0.272, 0.264, 0.267, 0.269, 0.289, against 0.225 to 0.226 for the chain at that e^{−f}: an
  unresolved horizon reads high, by 0.004 to 0.065, but not to their 0.296 to 0.42. So the test shows the direction and not the size, and it is not their scheme. The discrepancy with SH02 therefore **remains**.

## What critical_ghost.md §1.3 should now say (proposal for the coordinator)

Replace the last two paragraphs of "Is the threshold black-hole mass finite?" (from "The **0.30 a against 0.22 a discrepancy is unresolved**") by:

> *Added from #41* (`wormhole_evolution/settled_mass.md`; *numerical evidence*). All three papers define M as a Misner–Sharp mass with G = 1 in units of the throat radius, so no unit factor separates them;
> m_AH = r_AH/2 on the apparent horizon. In the symmetric family S at D = 4 the mass on the apparent horizon falls from 0.494 a at the first trapped cell and settles at 0.22308 a (r_AH = 0.4462 a), the same
> for p*ₕ − δ with δ = 1e-8, 1e-6 and up to 1e-3, from h = 0.02 to 0.0025 (second order), to e^{−f} = 1.6e6 on the horizon, 9.4 time units after formation. That is [GGS09b]'s r_AH = 0.45. [ST17]'s E_f
> equals the static throat's own energy (E_f/E_throat = 1.00 at n = 4, 5, 6), a horizon at r = a₀ read 2 to 3 time units after it forms: a formation value, which is why it does not depend on the injected energy.
> [SH02]'s 0.30 a is not reproduced with their own data (E₀ and x_H⁻ are, M is 0.223 to 0.26): it equals this code's M_AH about two time units after formation, and the mass is still falling then. The
> discrepancy is resolved for [GGS09b] and [ST17] and remains for [SH02]; a resolution limit of their reading is a conjecture and not shown.

and in §4 item 1 and the Summary, "0.30 a [SH02] against 0.22 a [GGS09b]" becomes "0.22308 a settled here, [GGS09b]'s 0.22; [SH02]'s 0.30 a not reproduced".

## Limits

- D = 4 and family S (plus SH02's data) only. D = 5, 6 are measured in #246 ([settled_mass_d56.md](settled_mass_d56.md)): S settles at 0.14587 a² and 0.10792 a³, G at 0.2640 a² and 0.8629 a³ (δ = 1e-8). The asymmetric static members that one-sided families may approach (#42, #225) are not run. Family G is measured in #253 ([settled_mass_g.md](settled_mass_g.md)): 0.2355 a at δ = 1e-8, 0.0125 a above S's, and it depends on δ.
- "Settled" means to e^{−f} = 1.57e6 (9.4 time units after formation). The mass is not followed to the event horizon or to I⁺.
- The AN control uses one h (0.02) and bar 3e-4, set above the observed 2.0e-4 after the development run. The RC6 band [1.5, 2.5] is wider than #39's, as the raw README explains.
- The tail of the SH02 explanation (their grid's resolution) is a conjecture. The uniform-grid test uses this code's kernel and not theirs.
- A log was overwritten by a launcher slip (the first `sh02` stdout at h = 0.01, `runs/sh02_h0.01.log`); its content is in `runs/sh02_h0.01.json` and the rerun log.

## Recommended tasks

1. *SH02's reading with a second code (#41).* What: write a second, independent double-null scheme (coordinator's), run SH02's Table I data on a uniform grid of 9601 points over x⁺ = [0, 20] without zoom, and read E(x⁺, x⁻)
   at x⁺ = 12, 16, 20 as their Fig. 7 does. Why: this decides whether 0.296 to 0.42 is an unresolved-horizon reading. Acceptance: M at the common point of the three curves within 5 % of SH02's six printed values, or
   a stated mismatch; a planted Tangherlini M recovered; a mutation of the grid size changes M. Depends on: this note.
2. *The settled mass beyond the limit.* What: K = 28 (e^{−f} to 4e8) and dxp = 18 at h = 0.005, δ = 1e-8, same checks. Why: the flux returned once between 4.7 and 7.5; RC8 only says it is below 1e-6 at 9.4. Acceptance: plateau mean unchanged to 1e-5 out to
   x⁺ − x⁺_f = 15; P1 to the new limit. Depends on: this note.
3. *Family G and the asymmetric static member (for #42, #225).* What: the same chain for G (threshold at a drifting throat). Why: whether the settled mass depends on the static member. Acceptance: M_AH(t) converges at three h; PL1 and PL2 pass unchanged; compare with 0.22308. Depends on: #225.
4. *Done in #246.* *D = 5, 6.* What: the chain at n = 3, 4 (M = r^{n−1}/2, a planted Tangherlini for each). Why: [ST17]'s n = 5, 6 rows. Acceptance: planted M recovered at second order. Depends on: this note.
5. *The event horizon.* What: integrate the outgoing null ray that stays on the horizon and compare r_EH/2 with M_AH. Why: GGS09b show the apparent horizon lies outside the event horizon; the settled mass is r_EH/2. Acceptance: r_EH agrees with r_AH to the plateau spread. Depends on: this note.

## Files

`settled_mass.py` (commands `pstar`, `pstar_an`, `job`, `planted`, `analyze`, `scan`, `sh02`, `sh02_uniform`), `settled_mass_output.txt`, this note; raw `negative-mass-raw/I41-rpi5`
(`runs/`: p*ₕ, the eight series, scan, sh02, sh02_uniform; `mut/`: the 15 mutated analyses; `dev/`: development logs and scripts; `papers.sha256`).

## Coordinator's review

Review of branch `wip/i41-threshold-masses` at dd4cad9 and raw `I41-rpi5` at 9d6082b, recorded in negative-mass-raw
`I41-mac-review` (criteria 8be196a; amendment A1; runs 020171c). On the Mac, nice 10.

- **Scope.** The branch adds `settled_mass.py`, its output and this note, nothing else.
- **Rerun.**
  - `analyze` gives 13 of 13, the same lines.
  - The job h = 0.01, δ = 1e-6 again agrees with the runner's M_AH(x⁺) to 1e-9 for the first half of the run. Deep in the
    zoom it agrees to 1.6e-5 at the end (0.223109 against 0.223125).
  - So the same code on another machine reproduces the settled value to about 2e-5. That is above the finest-h difference
    this note quotes (2.4e-6), and it is why the paper and `critical_ghost.md` quote 0.2231 a.
- **Independent check** (`ind41.py`). The review's own half-domain double-null march, sharing only the family's
  definition, reads the apparent horizon row by row on a uniform grid.
  - At e^{−f} = 1.335, 1.936, 4.074 and 11.96 it gives M_AH = 0.43726, 0.37797, 0.30315 and 0.25010 at h = 0.00125,
    within 1.2e-4 of this note's table. The horizon's x⁻ agrees to 4e-5.
  - Its first pair (h = 0.005, 0.0025) did not resolve the horizon at e^{−f} = 11.96 (K1 not met, amendment A1).
  - A uniform grid cannot reach the settled value (e^{−f} ≥ 1e5), so the plateau rests on this note's zoom chain and its
    planted Schwarzschild (to 1.7e-6).
- **Process.** Two choices were made after development runs and before the runs of record, without a ruling:
  - the declared limit moved from e^{−f} = 6144 to 1.57e6;
  - the AN bar was set at 3e-4, above an observed 2.0e-4.

  Both are recorded in the raw README. The first is justified by the flux that had not died. The second makes AN a weak
  test, and it is read as such.
- **Applied.**
  - `critical_ghost.md` §1.3 and §4 carry the proposed text, with 0.2231 a.
  - Paper III's threshold section states the settled mass. `paper_numbers.py` checks it (0 failures).
- **Follow-ups.** Tasks 1 and 4 are #245 and #246; tasks 2 and 5 together are #247. Task 3, family G, is added to #225.

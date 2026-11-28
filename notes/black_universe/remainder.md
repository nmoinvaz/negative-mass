# The remainder of smooth c = −1 data does not decay: the late growth at fixed ρ is c₂e^{κv}Y_∞, not c₂e^{κv}R (#322)

Code [remainder.py](remainder.py), figure script [remainder_graph.py](remainder_graph.py), figure [remainder.svg](remainder.svg), output of record [remainder_output.txt](remainder_output.txt), assembled from negative-mass-raw `I322-mac/runs/r1`; the criteria,
written before the runs of record, are in `I322-mac/README.md`. Builds on #302 ([horizon_basis.md](horizon_basis.md)), #262 ([cminus1_mode.md](cminus1_mode.md)) and #248 ([horizon_sector.md](horizon_sector.md)). Status words as in AGENTS.md: every number from the solver is
*numerical evidence*; the closed form of Sec. 1 is #262's (*checked symbolically* there, re-checked here to 1e-9). No sentence here says a result is new (the literature was not read for it).

## Answer

**1. Y_rem does not decay. It grows at the rate κ, with a profile that #262's closed form gives.** At fixed ρ (Π = −VU = −0.1, 0.5, 1, that is ρ = −0.3114 inside the horizon, 0.6819, 1.0000) the remainder of each of the 11 evolved data with c₂ ≠ 0 (the classes A0, A1R, A1L, L, G, Lm, Gb of #248 and #302, and 4 random data) has `Y_rem/R = aV + b + d/V` on V = e^{κv} ∈ [8, 24] (v from 3.3 to 5.0) with

```
a / c₂ = Q(ρ) = −J(ρ)/I₀ = Y_∞(ρ)/Y_ref(ρ) − 1 :   +0.29606 (Π = −0.1),   −0.57371 (Π = 0.5),   −0.73271 (Π = 1)
```

to 4.3e-4 at worst (class G, Π = 0.5; typically 1e-4), Richardson of h = 0.000625, 0.0003125, and the pair 0.00125, 0.000625 differs by ≤ 4.2e-5. Here `Y_∞ = Y_ref(1 − J/I₀)` is #262's solution of the equation at W = κ that is regular on the horizon and outgoing at infinity,
`J(ρ) = ∫₀^ρ 2Kr e^{−2H}/(AR²)dr`, I₀ = J(∞) = 1.0149036099284. The rate of Y_rem is ακ with the window exponent α of |Y_rem| = 0.940 (G) … 0.997 (Lm) at the exterior ρ; between the two Richardson pairs and the two finest raw h it changes by ≤ 6e-5
(exterior) and ≤ 4e-3 (interior ρ), so the rate is (0.94 … 1.00)κ over the window and tends to κ in the model (V-term plus a constant). It is not 0 and it is not negative.

**2. So the late growth at fixed ρ is `c₂e^{κv}Y_∞(ρ)`, not `c₂e^{κv}R(ρ)`.** The answer to #302's open question is no. Decomposition of smooth data at fixed ρ, for V ≤ 24:

```
Y(V, U)
├─ c₁ R(ρ)                     static mass mode (constant in v)
├─ c₂′ U R(ρ)                  decaying partner, ∝ 1/V
├─ c₂ V R(ρ)                   #262's mode, rate κ, profile R
└─ Y_rem
   ├─ −(c₂/I₀) V R(ρ) J(ρ)     #262's second solution Y_b = Y_ref J at W = κ: rate κ, vanishes on both horizons
   └─ y₁ f₀(ρ) + O(1/V)        V-independent over V ≤ 24; f₀: the W = 0 solution that vanishes at the horizon
```

Both Y_ref and Y_∞ have the horizon values c₁ + c₂V, so the horizon law and c₂ do not tell them apart; the profile at fixed ρ does. #262 says the condition at infinity picks Y_∞ in the two-dimensional space of solutions at W = κ; the evolution picks it too.

**3. The V-independent part** has the profile of the solution f₀ of #262's equation at W = 0 with f₀ = ρ + … at the horizon (the free coefficient y₁₁ of a VU term, D6 of #248): `y₁ = bR/f₀(ρ)` agrees at the three Π to 4.1e-2 (L), 2.1e-2 (G), 1.0e-2 or better for the others. For data with c₂ = 0 (class R, Rb, which have c₂ at the truncation seed,
1e-17 and 1e-12) there is no growth, α = −0.07 … +0.10 against 1 (the aV term is ≤ 5 % of b at V = 24), and the profile is f₀ again (y₁ spread 4.1e-2, 3.5e-3). Whether the static part decays later is not shown (Sec. 6).

(Figure: [remainder.svg](remainder.svg). (a) Y_rem/(c₂VR) of class G against V at the three ρ, with −J/I₀ dashed; (b) a/c₂ of the 11 data at the three ρ against the curve −J(ρ)/I₀.)

| acceptance (issue #322) | result |
|---|---|
| Y_rem's rate at three fixed ρ, with its error, at two h | rate = ακ, κ = 2/π. A0: α = 0.9912, 0.9870, 0.9868 at h = 0.000625 and 0.9881, 0.9870, 0.9868 at h = 0.0003125 (Π = −0.1, 0.5, 1); G: 0.9449, 0.9429, 0.9399 and 0.9447, 0.9429, 0.9399. All 13 data in Sec. 2 and `remainder_output.txt`. Error (largest of the two-pair change, the two-h change and the fit standard error): ≤ 6e-5 exterior, ≤ 4e-3 interior |
| A planted Y_rem with a known rate is recovered through the same fit | exact growing mode, 99 % of c₂ subtracted: α = 1.00000 (1.00000, 0.99983 interior), a/c₂ = 0.010000 (0.009993 interior); decaying mode left in: α = −1.00000, amplitude 1.000000 \|c₂′\|Π; ε Π V^α added at the nodes, α = 2, 1, 0, −1: α to 1e-4, amplitude to 1e-4 (bar 1e-3); the exact growing mode has Y_rem = 0 to 2e-9 (exterior), 2e-6 (interior) after Richardson, falling by 4.00 per halving of h |
| Whether the late growth is c₂e^{κv}R | no: it is c₂e^{κv}Y_∞; a/c₂ = −J/I₀ ≠ 0 to 4.3e-4 |
| each check can fail; mutations | 14 of 14 fired sets equal the declared ones (Sec. 4) |

## 1. The closed form (K1)

#262: with x = ∫dρ/A, e^{2κx} = ρe^{2H} (H(ρ) = ∫₀^ρ(K/A − 1/2r)dr analytic, H′(0) = π/8), R = √(1 + ρ²), the equation `Y″ = PY′ + QY` at W = κ has the solution Y_ref = Re^{κx} and, by Abel (Wronskian 2Kρ²/A), the second solution `Y_b = Y_ref J`,
`J(ρ) = ∫₀^ρ 2Kr e^{−2H}/(AR²)dr` (J = ρ + O(ρ³)). The combination with no e^{κx}-growth toward infinity is Y_∞ = Y_ref(1 − J/I₀). K1 re-derives J by DOP853 from |ρ| = 1e-4: I₀ = 1.01490360993 to 3.2e-13 of #262's value (K1a), residual of Y_ref and of Y_ref J
over the size of the terms 2.8e-16 and 1.7e-16 at ρ = −0.8 … 5 (K1b), 1 − J(20)/I₀ = 2e-12 (K1c). In Kruskal terms the fixed-ρ profile is e^{κt}Y(ρ) = V·R(1 + Q) with V = e^{κ(t+x)}.

## 2. The measurement

Smooth data: #302's classes A0, A1R, A1L, L, G, Lm, Rb, Gb, class R of #248 and 4 random data (seed 3221). Each is evolved in #248's Kruskal chart through S, V ∈ [−0.1, 24], U ∈ [−2, 0.06], at h = 0.00125, 0.000625, 0.0003125 by a streaming march (`march_stream`: `hs._cellK`'s cells, two columns in memory; the horizon line
agrees with `hs.march_K` to 3e-15, 8-node windows exactly; the interior corner is cut at Π = −0.3 because the interior ends at −0.37). At each V a window of 8 nodes around `U* = −Π/V` gives `Y/R − (c₁ + c₂V + c₂′U)` by 8-point Lagrange interpolation, with (c₁, c₂, c₂′) read at S by #302's slice functionals.
Richardson `(4w(h/2) − w(h))/3` pointwise on the common nodes. The fits are over V ∈ [8, 24]: `ln|Y_rem| = α ln V + const` (rate ακ) and `Y_rem/R = aV + b + d/V`. The error of α is the largest of the change between the two Richardson pairs, between the two finest raw h, and the standard error of the log fit.
**An error bar here is a spread over h and fit choices, not an enclosure.**

Table of record (Richardson of h = 0.000625, 0.0003125; `a/c₂ − Q` for the exterior Π = 0.5; Π = −0.1 and 1 in the output file):

| datum | c₂ | Q fit (Π = 0.5) | −J/I₀ | fit − pred | α (h = 0.000625 / 0.0003125) at Π = 0.5 | y₁ (Π = −0.1, 0.5, 1) |
|---|---|---|---|---|---|
| A0 | 2.837763e-02 | −0.57360 | −0.57371 | +1.1e-04 | 0.9870 / 0.9870 | −4.79e-3, −4.81e-3, −4.82e-3 |
| A1R | 1.337668e-03 | −0.57362 | −0.57371 | +9.0e-05 | 0.9911 / 0.9911 | −1.507e-4, −1.513e-4, −1.518e-4 |
| A1L | 1.541571e-02 | −0.57357 | −0.57371 | +1.4e-04 | 0.9826 / 0.9826 | −3.53e-3, −3.54e-3, −3.55e-3 |
| L | −9.942330e-04 | −0.57348 | −0.57371 | +2.3e-04 | 0.9558 / 0.9558 | 5.67e-4, 5.59e-4, 5.44e-4 |
| G | −9.942330e-04 | −0.57329 | −0.57371 | +4.3e-04 | 0.9429 / 0.9429 | 7.44e-4, 7.41e-4, 7.29e-4 |
| Lm | −2.725805e-03 | −0.57370 | −0.57371 | +8.2e-06 | 0.9969 / 0.9969 | 1.105e-4, 1.110e-4, 1.111e-4 |
| Gb | −2.725805e-03 | −0.57382 | −0.57371 | −1.1e-04 | 0.9907 / 0.9907 | 3.62e-4, 3.62e-4, 3.61e-4 |
| rand0 … rand3 | 7.9e-3, −1.04e-2, −1.85e-2, 1.19e-2 | −0.57355 … −0.57384 | −0.57371 | −1.2e-4 … +1.6e-4 | 0.9784 … 0.9919 | spread ≤ 1e-2 |
| Rb, R (c₂ = 0) | −1e-12, 6e-18 | (none) | | | 0.083 / 0.073; −0.053 / −0.053 | 2.51e-4, 2.51e-4, 2.50e-4; 1.78e-4, 1.81e-4, 1.85e-4 |

**An independent readout (D3).** The transverse slope of Y on the future horizon, `∂_UY(V, 0) = ∂_Uδr + (G₀V/√2)δφ`, from five nodes across U = 0, uses no fixed-Π window, no interpolation and none of the fits above. With Y/R = c₂V(1 + Q(ρ)) + y₁f₀ and ρ = −G₀VU it grows as `AV² + BV + C`
with A = c₂G₀/I₀, B = −y₁G₀: A/(c₂G₀/I₀) = 0.99918 … 1.00006 for the 11 data (bar 1e-2), B against −G₀ȳ₁ within 8 % (not a check). It shares with the fixed-Π readout the solver's output and G₀ = 2.47254596; it shares with the prediction nothing but G₀ and I₀.
**Shared between the prediction and the solver:** A(ρ), R(ρ), κ and the map ρ(Π) of #248's chart; not shared: the characteristic evolution, the interpolation and the fits.

## 3. Planted cases (P1 – P3)

The exact modes of #302 (β = 2e-3, c₂ = 1e-3 and c₂′ = −1e-3), as ray data of the same streaming march, V ∈ [0.05, 24] and U ∈ [−2, −0.05], h = 0.00125, 0.000625, 0.0003125, with c₁, c₂ (c₂′) read at the horizon crossings V₀ = 0.5, 1 (U₀ = −0.5, −1).

| id | planted | known answer | result |
|---|---|---|---|
| P1 | the growing mode | Y_rem = 0 | max \|Y_rem/R\|/(c₂V) on [8, 24]: Π = 0.5: 2.7e-4, 6.8e-5, 1.7e-5 (ratios 4.00, 4.00); Π = 1: 5.7e-5, 1.4e-5, 3.5e-6; Π = −0.1: 2.2e-2, 5.4e-3, 1.4e-3 (3.98, 4.00); Richardson 2.4e-9, 4.2e-10, 2.2e-6 |
| P2a | the growing mode, 99 % of c₂ subtracted | α = 1, a/c₂ = 0.01 | α = 1.00000, 1.00000, 0.99983; a/c₂ = 0.010000, 0.010000, 0.009993 (Π = 0.5, 1, −0.1) |
| P2b | the decaying mode, c₂′UR not subtracted | α = −1, amplitude \|c₂′\|Π | α = −1.00000 and amplitude ratio 1.000000 at Π = 0.5 (V ∈ [4, 9]), 1 (V ∈ [4, 18]) |
| P3 | ε Π V^α at the nodes of the exact mode, 2e-4 at V = 16, Π = 1 | α = 2, 1, 0, −1 | α within 5e-5 (printed to four decimals) and amplitude ratio 1.0000 at Π = 0.5, 1 |

The truncation floor of the exact mode falls by 4.00 per halving of h at every Π (second order); the interior floor is 80 times the exterior one (the interior curve at Π = −0.1 has U = 0.1/V ≈ 4e-3, close to the corner where the interior ends), which is why the interior α carries the larger error.

## 4. Checks and mutations

| id | what | result |
|---|---|---|
| K1a, K1b, K1c | the closed form of #262 | 3.2e-13; 2.8e-16, 1.7e-16; 2e-12 |
| P1, P2, P3 | planted zero and planted rates | pass (Sec. 3) |
| D1 | a/c₂ against −J/I₀ to 2e-3, pair change ≤ 1e-3, α ∈ [0.85, 1.02], 11 data | worst 4.3e-4, 4.2e-5, α 0.940 … 0.997 |
| D2 | y₁ = bR/f₀ equal at three Π to 0.10 | worst 4.1e-2 (L) |
| D3 | \|A/(c₂G₀/I₀) − 1\| ≤ 1e-2 | worst 8.2e-4 |
| D4 | c₂ = 0 data: α ≤ 0.3, y₁ spread ≤ 0.10 | α ≤ 0.095, spread ≤ 4.1e-2 |

All 10 pass in the runs of record (`derive` 3 of 3, 8 s; `check` 7 of 7, 9 s on the cache, 1.2 GB; `warm` 1837 s and 2045 s in two processes; code 50e1d53). 14 mutations, each with the checks it must fire declared before the runs of record (raw README; computed in development runs
on the cache of the 13 data; D3's two mutations on class G and A0 at the two coarser h): the fired sets equal the declared ones (14 of 14; `remainder_output.txt`).

| mutation | fires | effect in the run of record |
|---|---|---|
| k_W, k_J | K1b; K1a, K1b, K1c | W = 1.01κ: residual 6.5e-3; e^{+2H}: I₀ = 2e217 |
| no_c2, c2_sign | P1, P2, P3, D1 | Y_rem/(c₂V) = 1.000, 2.000; D1 1.0, 2.0 |
| no_c1 | D2 | y₁ of A0 −1.8e-2, 2.0e-3, 3.6e-4 |
| c2p_V | D1, D4 | D1 3e-2 … 9e-2; Rb α ≈ 0.94 … 1.02 |
| pi_off | P2, P3, D1 | amplitudes 1.050; D1 1.8e-2 |
| interp_near | P1, P2, D1, D2 | P1 Richardson 1.9e-4; D1 pair change 2.3e-2 |
| fit_semilog | P2, P3, D1 | P2a α = 0.066; D1 α = 0.065 |
| J_sign, plant_alpha | D1; P3 | 1.5; α + 0.1 |
| hp_sign, g_nophi | D3 | A/(c₂G₀/I₀) = 4.0 (A0, two coarser h in development), −1.51 (A0), 0.75 (G) |
| f0_W | D2, D4 | y₁ spread 0.5 |

Blind spots, kept: `no_c1` is not seen by P1 – P3 or D1 (the exact modes carry c₁ ≈ 1e-9; a constant does not change a/c₂); `no_c2` and `c2_sign` are not separated; `interp_near` is invisible to P3; `pi_off` leaves every exponent unchanged and is seen through the amplitudes only.

## 5. Proposed text (the files are not touched here)

| file, line | now | change |
|---|---|---|
| `horizon_basis.md` Sec. 3, last sentence of "What survives" | "That `Y_rem` decays, so that the late growth is c₂VR, was not tested." | *Tested in #322: it does not decay. At fixed ρ, Y_rem → −(c₂/I₀)·V·R J(ρ(Π)) + y₁f₀(ρ) (V ≤ 24), so the late growth is c₂VY_∞ with Y_∞ = Y_ref(1 − J/I₀).* |
| `horizon_basis.md` Sec. 7 first bullet; Sec. 8 task 1; the review's last bullet | "That `Y_rem` of smooth data decays …"; the task; "What is open: whether the remainder Y_rem of smooth data decays" | *Done in #322* (answer: no) |
| `horizon_sector.md` 12–13 (after #302's correction) | "… a growing solution of rate κ, regular on the whole slice, excited by regular data …" | append: *At fixed ρ the late growth is c₂e^{κv}Y_∞(ρ) with Y_∞ = Y_ref(1 − J/I₀) of #262, not c₂e^{κv}R (#322).* |
| `horizon_sector.md` "What this does not show", bullet "The late-time growth …" | "The rate κ is #232's; here only the horizon law and its coefficient are tested." | *Done in #322 to V = 24 (v = 5.0): rate κ at fixed ρ, profile Y_∞.* |
| `cminus1_mode.md` Sec. 3 ("What this means"), after "the condition at infinity picks Y_∞ in it" | | *#322: evolved smooth data reach Y_∞ at fixed ρ: the growth coefficient a/c₂ of Y − (c₁ + c₂V + c₂′U)R equals Y_∞/Y_ref − 1 = −J/I₀ to 4.3e-4 (11 data, V ∈ [8, 24]).* |
| `paper/4-black-universes/main.tex` 471–475 ("The decaying partner"), last clause | "a remainder that vanishes on both horizons; whether the remainder decays is open (\texttt{notes/black\_universe/horizon\_basis.md})" | "a remainder that vanishes on both horizons, which does not decay (below)" |

**Paper IV**, after the paragraph "The decaying partner" (the paper's statement at line 462 that the family linearizes to the solution that grows toward infinity stands: that solution is Y_ref):

```
\paragraph{The remainder (\#322).} The remainder does not decay. At fixed $\rho$ it grows at $\kappa$ with the coefficient $a/c_2=-J(\rho)/I_0$, $J=\int_0^\rho 2Kr\,e^{-2H}/(AR^2)\,\dd r$ and $I_0=J(\infty)=1.01490361$,
which is the profile of $Y_\infty=Y_{\rm ref}(1-J/I_0)$, the solution regular on the horizon and outgoing at infinity. Evolved smooth data (eleven, $V=e^{\kappa v}$ to 24) give $a/c_2$ within $4.3\times10^{-4}$ of it at
$\rho=-0.311$, $0.682$, $1$ (Richardson, $h=6.25\times10^{-4},3.125\times10^{-4}$). The late growth is $c_2e^{\kappa v}Y_\infty$, not $c_2e^{\kappa v}R$: both have the horizon values $c_1+c_2V$ and differ by the second solution
$Y_b=Y_{\rm ref}J$ of the same equation. What is left is independent of $V$ over this range and has the profile of the static solution that vanishes on the horizon.
```

## 6. What this does not show, and what disagrees

- **V > 24 (v > 5.0).** The fits are over v from 3.3 to 5.0. The model aV + b + d/V has a residual of 1e-6 to 5e-6 of the signal; whether the V-independent part y₁f₀ decays later is not shown. #232 (Paper IV, Sec. on the horizon sector) found that its combination of three pulses with c₁ = c₂ = 0 decays at the rate −0.1768 at x = 0, t = 14;
  class R and Rb here (c₁ = c₂ = 0) are V-independent to 5 % over v ≤ 5.0 at ρ ≥ 0.68. The two are not in conflict (t = 14 is V ≈ 7·10³; the ρ differ), but #232's decay is not reproduced here either. #232's chart reaches late times; this one does not.
- **Why the evolution picks Y_∞.** #262 says the condition at infinity picks Y_∞. Here the data carry no incoming flux of their own on the ray U = −2 beyond what the classes specify (L, G, Lm: none), yet A0, A1R, A1L, whose ray data do not vanish, give the same Q; no derivation is given.
- **Planted Y_b.** The solver was not given an exact Y_b mode (no closed form of its fields); the planted rates are those of Y_ref's time reverse and of known functions at the nodes.
- **Interior.** The interior curve Π = −0.1 has a truncation floor 80 times the exterior (Sec. 3); its a/c₂ agrees to 2.5e-4 after Richardson, its α carries 4e-3.
- **Other members.** At c ≠ −1 the horizon is not at the area minimum (ρ_h ≠ 0, #301) and the single growing mode is a regular quasinormal-type solution, so f₀ and Y_∞ do not carry over.
- **Nonlinear.** Nothing about #217's lifetime law, whose growth is that of the exact family (Y_ref).
- Literature: nothing was read.

## 7. Recommended tasks

1. **The static part at late times.** *What:* evolve class R and Rb, whose Y_rem is only y₁f₀, to V ≈ 10³ (v ≈ 11) in a chart that resolves it (the (p, σ) chart of #232 with #302's data, or a hyperboloidal Kruskal grid) and fit the late rate of y₁f₀. *Why:* #232's c₁ = c₂ = 0 combination decays at −0.1768 and R, Rb do not to v = 5.
   *Acceptance:* the rate of y₁f₀ at three ρ with its error at two h, planted (an exact decaying f₀-type solution or a known tail) recovered, and a mutation that drops the f₀ subtraction fails. *Depends on:* #302, #232.
2. **y₁ from the data at S.** *What:* the free coefficient y₁₁ of D6 (the Π-linear term of Y at S) from the slice data, and its relation to y₁ = bR/f₀. *Why:* y₁ is a number of the data that the three of #302 do not read. *Acceptance:* y₁ from the slice equals y₁ of Sec. 2 to 1e-3 on the 11 data; a wrong f₀ fails. *Depends on:* #302 task 2, #248 D6.
3. **A planted Y_b.** *What:* exact double-null ray data for V Y_b (and for V Y_∞) from the ODE solution and a gauge vector as in #302, evolved through the solver. *Why:* a planted remainder of rate κ that is a solution of the system, not a function at the nodes. *Acceptance:* the fit returns a/c₂ = −1/I₀·(…) for Y_b and Q(ρ) for Y_∞ to 1e-3; the mutation `J_sign` fails it. *Depends on:* #302, #262.

## 8. Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, two processes for `warm`:

```
remainder.py derive [--mutate=k_W|k_J]                              # K1, 8 s
remainder.py warm --cache=DIR [--items=A0,G,...,exact_g,exact_d]    # evolves the items at three h, about 280 s per item, DIR 0.9 GB
remainder.py check --cache=DIR [--only=P1,..,D4] [--mutate=no_c2|c2_sign|no_c1|c2p_V|pi_off|interp_near|fit_semilog|J_sign|plant_alpha|hp_sign|g_nophi|f0_W]   # 7 checks, 9 s on the cache
remainder.py table --cache=DIR                                      # the numbers of this note
remainder.py dump OUTDIR --cache=DIR                                # the Richardson series of every datum
remainder.py mutations DIR PREFIX                                   # fired sets of the stored logs DIR/PREFIXname.log against the declared ones
remainder_graph.py SERIES_DIR                                         # remainder.svg
```

`remainder_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I322-mac/runs/r1`, code commit 50e1d53), in that order. The cache of evolved readouts is not committed (0.9 GB); `warm` regenerates it.

# Radial stability of the phantom black universes

November 2026, #212 (branch `wip/i212-black-universe`; raw `I212-rpi5` in negative-mass-raw). Code `stability.py`, output
`stability_output.txt`, graph `growth_rate.svg` (`rate_graph.py`). Literature note: #211. Context: `../critical_ghost.md` §1.1.

Status words as in AGENTS.md. Every count below is *numerical evidence* (a converged computation); the algebra is *checked symbolically*
(sympy, inside `stability.py`). Nothing here is a proof.

## Result

- Nineteen of the twenty black universes tested (c from −10 to −0.01 in units b = 1, c ≠ −1) have **exactly one growing radial mode**; the
  member c = −1, whose horizon sits at the area minimum, has **none**. **Update from #232:** none under the condition used here, u = Y/ρ → 0 at the horizon, which sets Y = 0 on it. For perturbations with Y regular and nonzero on the horizon, the linear double-null response grows at κ (Y_h ∝ e^{κv}; numerical evidence, and the horizon exponents 1 ± W/A₁ checked symbolically in #232's review). Whether such perturbations are free in physical Cauchy data is #248. R1 and R2 agree to 1.2e-8 or better at every member, and an
  independent spectral solver gives the same counts. **Update from #248** ([horizon_sector.md](horizon_sector.md)): such perturbations are physical. The coefficient c₂ = ∂_V δr of the growing sector on the horizon is a free number of regular Cauchy data on a slice through the horizon; the solution Y = c₂ V G(Π) is smooth across the bifurcation surface and into the interior; 20 of 20 random regular data have c₂ ≠ 0. So the linear problem at c = −1 has one growing direction of rate κ, like the other 19 members, and #212's count of none is the count for c₁ = c₂ = 0 (u = Y/ρ → 0), which data supported away from the horizon crossing satisfy (checked symbolically and numerical evidence; the coordinator's transport along the data ray reproduces c₂ of a pulse to 2.3e-8). The mechanism is general ([horizon_lemma.md](horizon_lemma.md), #261, proved): on any non-extremal Killing horizon at an extremum of the areal radius, δr = c₁ + c₂V along the horizon with c₂ gauge invariant, and such a horizon needs A not a function of the areal radius alone.
- This is the verdict of [BKZ12] member by member: no disagreement at the 20 members, eight of which they test in the time domain.
- The stable member is isolated. For c → −1 from either side the growth rate tends to the surface gravity κ = A′(ρ_h)/2 = 2/π ≈ 0.6366
  (b = 1), not to zero; at c = −1 the mode is gone. The limit is an exact bound state (Sec. 5).

## 1. The family

[BF06] (gr-qc/0511109, PRL 96, 251101) and [BKZ12] (arXiv:1205.2224, PRD 86, 024028) use a phantom scalar (ε = −1) with the potential that
supports `ds² = A dt² − dρ²/A − R² dΩ²`:

```math
R = \sqrt{\rho^2 + b^2},\qquad \varphi = \sqrt2\arctan(\rho/b),\qquad
A = R^2\Big[\frac{c}{b^2} + \frac1{R^2} + \frac{3m}{b^3}\Big(\frac{b\rho}{R^2} + \arctan\frac\rho b\Big)\Big]
```

and `V(φ)` in [BF06] (14), [BKZ12] (4). Asymptotic flatness at ρ → +∞ is `2bc = −3πm` ([BF06] (13), [BKZ12] (5)). Units b = 1, so
`3m = −2c/π` and c is the only parameter.

| c | member |
|---|---|
| 0 | Ellis wormhole (V = 0, A = 1): the planted case |
| −1 < c < 0 | black universe, horizon at ρ_h < 0, the throat ρ = 0 inside the static region ρ > ρ_h |
| −1 | horizon at the throat, ρ_h = 0 (3πm = 2b) |
| c < −1 | black universe, ρ_h > 0, no throat in the static region |
| c > 0 | M-AdS wormhole; not computed here (**Update from #218:** one growing radial mode at each of 12 values, `polar.md`) |

Checked: the closed forms satisfy [BF06] (3)–(6) with residual 0 in sympy and ≤ 1e-38 at 40 digits (S1, S2). [BF06] (8) prints
`A′(h) = (ρ₀ − h)/R²(h)`; (7) gives `2(ρ₀ − h)/R²(h)`, which S3 confirms at 40 digits (so κ = (ρ₀ − h)/R²).

## 2. The perturbation equation

Derived in sympy from `R_μν − ε ∂_μφ ∂_νφ − V g_μν = 0` and the scalar equation, linearized about the generic static background, e^{Wt}
time dependence. Gauge **δφ = 0**: regular everywhere, since φ′ ≠ 0. [BKZ12] use δR = 0 (35), singular at R′ = 0.

```
tρ equation       δλ = (Y′ − A′Y/2A)/R′                         (Y = δR)
scalar equation   fixes δν′
θθ equation       master equation; the tt and ρρ equations then hold identically (residual 0)
```

For the family, with `dx = dρ/A` and `u = Y/ρ`:

```math
Y'' = \Big(\frac2\rho - \frac{A'}A\Big)Y' + \Big(\frac1{R^4} - \frac{A'}{A\rho R^2} - \frac2{AR^2} + \frac{W^2}{A^2}\Big)Y,\qquad
u_{xx} = (\hat V + W^2)u,\quad \hat V = \frac{2A^2}{\rho^2} - \Big(1+\frac1{R^2}\Big)\frac{AA'}{\rho} + \frac{A^2}{R^4} - \frac{2A}{R^2}
```

- At c = 0 it is `notes/wormhole_largeD/master.py`'s equation for D = 4 and reproduces W²a² = 1.39705243433 (Sec. 3).
- The throat is a regular singular point of Y with exponents 0 and 3, the wall 2A²/ρ² of V̂ (u ~ ρ⁻¹, ρ²). Both solutions are analytic there.

**Comparison with [BKZ12].** The gauge map `F = −φ₀′Y/R′` takes this equation onto the δR = 0 one (S4). Their potential (40), (41)
`V_eff = U + R_{r*r*}/R`, `U/A = ε φ′²/R′² (R²V − 1) + 2φ′R V_φ/R′ + ε V_φφ` equals V̂: symbolically `A·Ũ = rhs` of (41) (S4, ε = ±1) and numerically
to 1e-29 (S8), with `Ψ = −√2 u`. **The printed (36), (37), (38) do not agree with the derivation; (40), (41) do.**

| | printed | derived |
|---|---|---|
| (36), coefficient of the time derivative | `e^{−2μ₂−2ν}` = 1 | `1/A²` |
| (36), coefficient of δφ′ | `ν′ − 2μ₃′ + μ₂′` = −2R′/R | `ν′ − μ₂′ + 2μ₃′` = A′/A + 2R′/R |
| (38) | `δφ = Ψ e^{μ₃}` | `δφ = Ψ/R` |
| (37), first term | `ε(e^{2μ₃} − V)φ′²/μ₃′²` | `ε(V − e^{−2μ₃})φ′²/μ₃′²` |
| (37), V_φ term | `−2φ′V_φ/μ₃′` | `+2φ′V_φ/μ₃′` |

Reversing the signs of μ₂ and μ₃ in (36), (38) and in (37)'s V_φ term gives the derived forms, but not (37)'s first term. Not corrected;
nothing here uses them.

## 3. Boundary conditions and method

[BKZ12] (43), (44): quasinormal conditions `Ψ ∝ e^{±iωr*}`, r* → ±∞ (Ψ ∝ e^{−iωt}); ingoing at the horizon, outgoing at infinity. A growing
mode has ω = iW, so the solution decays at both ends and W² is an eigenvalue of a Schrödinger operator. Here:

- horizon: ingoing, `Y ~ (ρ − ρ_h)^{W/A′(ρ_h)}` (at c = −1, `ρ^{1+W/A′}`). Integration outward from ln(ρ − ρ_h) = ln z₀; the ingoing solution
  dominates outward, so the start vector cannot select the other one (a control that tried it was silent, see raw README);
- infinity: `u_x/u = −W − V̂/2W`, integrated inward;
- throat: no condition. Analytic continuation from −d to +d along the upper half circle in the complex ρ-plane; the monodromy test
  (relative imaginary part of Y afterwards) is at most 6.6e-12 over all throat members and both resolutions, so there is no logarithm.

A growing mode is a sign change of the Wronskian mismatch of the two solutions on W ∈ [0.01, 8] (units 1/b), 60 (R1) and 150 (R2) log-spaced
points, refined by Brent. The bound on the range searched is that interval. Each root's eigenfunction is integrated: no node at any member,
so each is the ground state (no root above it).

**Independent check (S6).** Chebyshev collocation of the same equation as a matrix problem (u-form for c ≤ −1, Y-form above), `N = 401, X = 45`
and `N = 801, X = 90`. It shares the background and the coefficients with the shooting code and nothing of the boundary conditions, the
throat crossing or the root finding, and finds complex eigenvalues too. Same counts at all 20 members (40 of 40), no complex eigenvalue with
Re λ > 1e-3. Rates against shooting: c ≤ −1, S-B 4.6e-7 at worst (c = −10, a long decay length) and 1e-9 or better elsewhere; throat members
(c > −1), 4.5e-3 at worst at both resolutions, since the Y-form near the throat does not converge with N.

**Planted (P1).** Through the same code at c = 0: one growing mode, `W²a² = 1.397052434006988` (R1) and `1.397052434332310` (R2) against
`1.39705243433350684` of `../wormhole_largeD/w2_output.txt`: errors 3.3e-10 and 1.2e-12, declared tolerances 1e-7 and 1e-9. With the kinetic
sign flipped (`--mutate=eps`), the same code gives three growing modes at c = 0.

## 4. The counts

W is in units 1/b. The bracket column is W in [BKZ12]'s units (m = 1/3), i.e. times 3m.

| c | ρ_h | κ | modes R1, R2 | W (1/b) | W (m = 1/3 units) | W/κ | nodes | R1 v R2 |
|---|---|---|---|---|---|---|---|---|
| −10 | 4.197082 | 0.116522 | 1, 1 | 0.0895350811 | 0.569998 | 0.76840 | 0 | 1.2e-10 ok |
| −3 | 1.118725 | 0.351374 | 1, 1 | 0.2827745372 | 0.540060 | 0.80477 | 0 | 3.6e-10 ok |
| −2 | 0.621184 | 0.470502 | 1, 1 | 0.3984523430 | 0.507325 | 0.84687 | 0 | 1.3e-09 ok |
| −1.5 | 0.339617 | 0.551682 | 1, 1 | 0.4941711675 | 0.471899 | 0.89575 | 0 | 1.2e-10 ok |
| −1.1 | 0.075763 | 0.620954 | 1, 1 | 0.6031905199 | 0.422403 | 0.97139 | 0 | 2.6e-10 ok |
| −1.01 | 0.007824 | 0.635123 | 1, 1 | 0.6331476235 | 0.407105 | 0.99689 | 0 | 9.3e-11 ok |
| −1.001 | 0.000785 | 0.636471 | 1, 1 | 0.6362712254 | 0.405468 | 0.99969 | 0 | 1.1e-10 ok |
| −1.0001 | 0.000079 | 0.636605 | 1, 1 | 0.6365849043 | 0.405303 | 0.99997 | 0 | 2.3e-10 ok |
| −1 | 0.000000 | 0.636620 | 0, 0 | none | none | – | – | 0.0e+00 ok |
| −0.9999 | −0.000079 | 0.636635 | 1, 1 | 0.6366546434 | 0.405266 | 1.00003 | 0 | 1.2e-08 ok |
| −0.999 | −0.000786 | 0.636768 | 1, 1 | 0.6369686164 | 0.405101 | 1.00031 | 0 | 3.0e-10 ok |
| −0.99 | −0.007884 | 0.638098 | 1, 1 | 0.6401216241 | 0.403439 | 1.00317 | 0 | 5.4e-09 ok |
| −0.95 | −0.040055 | 0.643811 | 1, 1 | 0.6544315613 | 0.395793 | 1.01650 | 0 | 7.9e-10 ok |
| −0.9 | −0.081826 | 0.650428 | 1, 1 | 0.6730206207 | 0.385612 | 1.03473 | 0 | 2.8e-09 ok |
| −0.7 | −0.271920 | 0.668151 | 1, 1 | 0.7558707523 | 0.336842 | 1.13129 | 0 | 2.5e-09 ok |
| −0.5 | −0.524073 | 0.660873 | 1, 1 | 0.8541005253 | 0.271869 | 1.29238 | 0 | 2.4e-10 ok |
| −0.3 | −0.926754 | 0.601300 | 1, 1 | 0.9700800261 | 0.185272 | 1.61331 | 0 | 1.8e-10 ok |
| −0.2 | −1.284188 | 0.532819 | 1, 1 | 1.0354140714 | 0.131833 | 1.94327 | 0 | 7.9e-11 ok |
| −0.1 | −2.024885 | 0.409506 | 1, 1 | 1.1059696680 | 0.070408 | 2.70074 | 0 | 8.7e-11 ok |
| −0.01 | −7.002156 | 0.140086 | 1, 1 | 1.1741182764 | 0.007475 | 8.38142 | 0 | 1.9e-10 ok |

- "R1 v R2" is the largest relative difference of the rates; the criterion was 1e-6.
- W falls monotonically with −c, from the Ellis value 1.18197 at c = 0 (1.1741 at c = −0.01) to 0.0895 at c = −10.
- In [BKZ12]'s units the rate tends to zero as c → 0⁻ and does not fall toward c = −1 (0.3856, 0.3958, 0.4034 at c = −0.9, −0.95, −0.99, then
  0.4053 at −1.0001). They describe both in words for Figs 7 and 8 and give no rates.
- The graph: `growth_rate.svg`.

## 5. Why c = −1 is alone

**Update from #248:** c = −1 is alone only in this section's sense, the L² count with u = Y/ρ → 0. With Y regular on the horizon it has a growing
direction of rate κ, whose coefficient is free in regular Cauchy data ([horizon_sector.md](horizon_sector.md)). What follows explains why the bound state of
the count disappears at c = −1, not why the member is stable.

**Update from #262** ([cminus1_mode.md](cminus1_mode.md)): at c = −1 the growing mode is Y = c₂ V R(ρ(Π)), that is
Y_ref = R e^{κx} with x the tortoise coordinate. It is an exact solution for every member at W = K = −2c/π (checked symbolically, and by the
coordinator's derivation from the linearized field equations in the static gauge, `I262-mac-review/ind262c.py`). The solution regular on the
horizon and decaying at infinity is Y_∞ = Y_ref − Y_b/I₀, I₀ = 1.01490361. The rate is smooth through c = −1:
W(c) = K(c) + (c + 1)/I₀ + O((c + 1)²), d(W/κ)/dc = 0.31402902 (numerical evidence; the third-order fit of the rates gives 0.31402903). The two-sided
0.31 of the table is this slope.

**Update from #269** ([exact_mode.md](exact_mode.md)): Y_ref is the linearization about the member of an exact nonlinear
solution: a time-dependent throat radius b(v)² = L + βe^{Kv}, with the scalar unchanged and F = A − 2ρλ′/λ + 1/λ² − 1 in the ingoing null chart
(checked symbolically for every member, and by the coordinator's brute-force Ricci check, `I269-mac-review/ind269.py`, 9e-40). For β < 0 it
collapses at finite advanced time with a diverging Ricci scalar; for β > 0 it expands. It is not asymptotically flat, and it linearizes to Y_ref,
not to Y_∞, so #250 stays open.

Near the horizon, with the scaling ρ = ρ_h t at fixed t and ρ_h → 0, V̂ in x tends to a fixed well (sympy, S7, S7c; the members c = −1.0001 and
−0.9999 follow it to 1.5e-3 and 1.3e-3, S7b, S7d):

| side | limit of V̂ | mode | W |
|---|---|---|---|
| c → −1 from below (ρ_h > 0) | `−(A₁²/2)/cosh²(A₁x/2)`, the ℓ = 1 Pöschl–Teller well | `u = 1/cosh(A₁x/2)` | κ = A₁/2 |
| c → −1 from above (ρ_h < 0) | `(A₁²/2)/sinh²(A₁x/2)`, with the throat's wall | `u = 1/sinh(A₁x/2)`, `Y = ρu` regular | κ = A₁/2 |

with A₁ = A′(ρ_h). The well stays at a fixed shape while it moves to x → −∞ like (1/A₁) ln|ρ_h|, so the bound state keeps W = κ. At c = −1 exactly
there is no well at finite x: V̂ → 0 smoothly at the horizon (the wall and the horizon cancel), and the computation finds no bound state (no root on [0.01, 8], no positive eigenvalue). The
rates agree: W/κ = 0.99969, 0.99997 at c = −1.001, −1.0001 and 1.00003, 1.00031 at c = −0.9999, −0.999. This is why the count is 0 at one point
and 1 at points arbitrarily close to it; it is not a numerical threshold effect.

## 6. Comparison with [BKZ12], member by member

[BKZ12] states (Conclusions 2): stable only at c = −1, unstable at c ≠ −1. For the members they show:

| c | modes found | [BKZ12] | their evidence | agree |
|---|---|---|---|---|
| −10 | 1 | unstable | Conclusion 2 | yes |
| −3 | 1 | unstable | potential, Fig 6 | yes |
| −2 | 1 | unstable | potential, Fig 6 | yes |
| −1.5 | 1 | unstable | potential, Fig 6 | yes |
| −1.1 | 1 | unstable | Conclusion 2 | yes |
| −1.01 | 1 | unstable | potential, Fig 5 | yes |
| −1.001 | 1 | unstable | time domain, Fig 5 | yes |
| −1.0001 | 1 | unstable | Conclusion 2 | yes |
| −1 | 0 | stable | time domain, Fig 5 (decays) | yes |
| −0.9999 | 1 | unstable | Conclusion 2 | yes |
| −0.999 | 1 | unstable | Conclusion 2 | yes |
| −0.99 | 1 | unstable | time domain, Fig 8 | yes |
| −0.95 | 1 | unstable | time domain, Fig 8 | yes |
| −0.9 | 1 | unstable | time domain, Fig 8 | yes |
| −0.7 | 1 | unstable | Conclusion 2 | yes |
| −0.5 | 1 | unstable | Conclusion 2 | yes |
| −0.3 | 1 | unstable | time domain, Fig 7 | yes |
| −0.2 | 1 | unstable | time domain, Fig 7 | yes |
| −0.1 | 1 | unstable | time domain, Fig 7 | yes |
| −0.01 | 1 | unstable | Conclusion 2 | yes |

No member disagrees. For the others the verdict is that of their Conclusion 2. They give no count of growing modes and no rates. What they
do show matches the computation in kind:
- their Fig 5 potential at c = −1 has a negative gap and the time-domain run decays (stable); ours finds no root and the spectral solver no bound state;
- they remark that approaching c = −1 the growth rate does not fall as one would expect near a stable member; here it tends to κ.

## 7. Members with exactly one growing mode

The 19 members with c ≠ −1, with the rate of the one growing mode:

- c = −10: W = 0.0895350811 per b (0.569998 in m = 1/3 units), e-folding time 11.169 b
- c = −3: W = 0.2827745372 per b (0.540060 in m = 1/3 units), e-folding time 3.536 b
- c = −2: W = 0.3984523430 per b (0.507325 in m = 1/3 units), e-folding time 2.510 b
- c = −1.5: W = 0.4941711675 per b (0.471899 in m = 1/3 units), e-folding time 2.024 b
- c = −1.1: W = 0.6031905199 per b (0.422403 in m = 1/3 units), e-folding time 1.658 b
- c = −1.01: W = 0.6331476235 per b (0.407105 in m = 1/3 units), e-folding time 1.579 b
- c = −1.001: W = 0.6362712254 per b (0.405468 in m = 1/3 units), e-folding time 1.572 b
- c = −1.0001: W = 0.6365849043 per b (0.405303 in m = 1/3 units), e-folding time 1.571 b
- c = −0.9999: W = 0.6366546434 per b (0.405266 in m = 1/3 units), e-folding time 1.571 b
- c = −0.999: W = 0.6369686164 per b (0.405101 in m = 1/3 units), e-folding time 1.570 b
- c = −0.99: W = 0.6401216241 per b (0.403439 in m = 1/3 units), e-folding time 1.562 b
- c = −0.95: W = 0.6544315613 per b (0.395793 in m = 1/3 units), e-folding time 1.528 b
- c = −0.9: W = 0.6730206207 per b (0.385612 in m = 1/3 units), e-folding time 1.486 b
- c = −0.7: W = 0.7558707523 per b (0.336842 in m = 1/3 units), e-folding time 1.323 b
- c = −0.5: W = 0.8541005253 per b (0.271869 in m = 1/3 units), e-folding time 1.171 b
- c = −0.3: W = 0.9700800261 per b (0.185272 in m = 1/3 units), e-folding time 1.031 b
- c = −0.2: W = 1.0354140714 per b (0.131833 in m = 1/3 units), e-folding time 0.966 b
- c = −0.1: W = 1.1059696680 per b (0.070408 in m = 1/3 units), e-folding time 0.904 b
- c = −0.01: W = 1.1741182764 per b (0.007475 in m = 1/3 units), e-folding time 0.852 b

Each mode is the ground state of its operator (no nodes), so each could serve as the critical
solution of a threshold of the kind the issue calls "type III" after [CHM99] (a critical solution with one unstable mode). Nothing is claimed here beyond that: whether
a threshold in phantom collapse approaches any of them is not tested.

## 8. `critical_ghost.md` §1.1, the row on [BKZ12]

The row says: "Phantom with a specific potential (M-AdS wormholes, 'black universes') | Unstable, except black universes whose horizon sits at
the area minimum". [BKZ12]'s abstract: "all configurations under study are unstable under spherically symmetric perturbations, except for a
special class of black universes where the event horizon coincides with the minimum of the area function". **The quote is fair.** Two points
need adding for it to be read correctly:
- the "class" is one member of a one-parameter family (c = −1, 3πm = 2b), not an open set;
- it is a knife edge: members at c = −1 ± 10⁻⁴ have a growing mode of rate κ to 3e-5, not a small one.
  **Update from #232:** the rate is κ at c = −1 as well, for perturbations with Y_h ≠ 0; the member is stable only against perturbations with Y_h = 0, such as waves incident on an unperturbed horizon (`linear.md`). **Update from #248:** those perturbations are free in regular Cauchy data, so the member decays only for data with c₂ = 0, a single condition at the horizon crossing, which includes data supported away from the horizon.

[BKZ12] show the instability of the M-AdS wormholes by regularized potentials (Fig 1) and time-domain runs (Fig 2); this work did not recompute them.

## 9. What `critical_ghost.md` and #211 should now say

Proposed wording, for the coordinator (files not touched here).

`critical_ghost.md` §1.1, the [BKZ12] row, result column: *Black universes (a one-parameter family with the throat radius b as the unit):
exactly one growing radial mode at each of 19 members tested, c = −10 … −0.01, rates 0.09 … 1.17 per b, none at the single member c = −1
whose horizon sits at the area minimum [BKZ12]; as c → −1 the rate tends to the surface gravity κ ≠ 0, not to zero (#212, numerical
evidence). M-AdS wormholes: unstable [BKZ12], not recomputed.*

#211, entry on [BKZ12] and on [BF06]'s suggestion: *The stability test of [BF06]'s solutions was computed (#212). Under radial perturbations
the black universes are linearly unstable, with one growing mode, except at one isolated member. The e-folding time 1/W runs from 0.85 b (c → 0⁻, the Ellis value) to 11 b (c = −10), with b the throat radius.
Nothing is shown about non-radial modes, the Kantowski–Sachs region, or formation from collapse.*

## 10. Checks and controls

| check | what | result |
|---|---|---|
| S1, S2 | background (3)–(6) on the closed forms, sympy and 40 digits | residual 0 and ≤ 1e-38 |
| S3 | A′(h) = 2(ρ₀ − h)/R² at the horizon | 5.5e-39 |
| S4 | linearization: tt, ρρ after the master equation, δR = 0 gauge, [BKZ12] (41), gauge map | all 0 |
| S5 | family coefficients against the generic equation | 3.0e-31 |
| S6 | spectral count and values | 40 of 40 pass |
| S7, S7b, S7c, S7d | the c → −1 limits (sympy) and the members | residual 0; 1.5e-3, 1.3e-3 |
| S8 | V̂ against [BKZ12] (40), (41), 30 digits | 1.1e-29 |
| P1 | planted Ellis | 3.3e-10 (R1), 1.2e-12 (R2) |
| R, V | R1 v R2, verdict against [BKZ12] | 20 of 20, 20 of 20 |

Fourteen `--mutate=` controls (`stability.py controls`), each fires on its target in the run of record: `bg_V`, `bg_A`, `bg_kappa`, `master_p`,
`bkz41_sign`, `family_q`, `veff_sign`, `pt_ell` (checks), `eps` (3 growing modes instead of 1), `throat_wall` (0 instead of 1), `rmax` (4.7e-7 against 1e-7),
`hz_start` (1.9e-4 against 1e-6), `scan_range` (0 modes), `spec_sign`. A control that was silent in development (the outgoing horizon exponent) was dropped and is
recorded in the raw README. The `throat_wall` control is the claim of [BKZ12] that the initial singular potentials are positive and that the
infinite wall hides the instability: with independent sides and a Dirichlet condition at the throat the Ellis wormhole has no growing mode.

## 11. What this does not show

- Only radial (ℓ = 0) perturbations in the static region ρ > ρ_h. The Kantowski–Sachs region beyond the horizon, polar ℓ ≥ 1 modes
  ([BKZ12] treat polar ℓ = 0 and the axial modes of the black universe at c = −1, Fig 9, Table I) and the M-AdS wormholes are not computed.
- Members of this one family and one potential only. Counts are for the 20 values of c listed; between them nothing is computed.
- The count rests on a sign-change search on [0.01, 8] with the node and spectral guards; no argument-principle count, no rigorous bound.

The development history, the controls that were silent or replaced, and the amendments are in the raw README (`negative-mass-raw/I212-rpi5/README.md`).

## Reproduce

From this directory, with `~/Source/negative-mass-venv314/bin/python`:

```
stability.py check; stability.py planted; stability.py count <c> --res=R1|R2 --out=member_<c>_<res>.json
stability.py report DIR; stability.py spectral DIR; stability.py controls; rate_graph.py
```

`stability_output.txt` is these outputs of record, assembled from negative-mass-raw `I212-rpi5`. About 14 minutes for the 40 count runs on
three processes, 27 s for `check`, 6 minutes for `controls`, 126 MB at most.

## References

- [BF06] K. A. Bronnikov, J. C. Fabris, "Regular phantom black holes", PRL 96, 251101 (2006), arXiv:gr-qc/0511109. Read in full (PDF and TeX).
- [BKZ12] K. A. Bronnikov, R. A. Konoplya, A. Zhidenko, "Instabilities of wormholes and regular black holes supported by a phantom scalar field",
  PRD 86, 024028 (2012), arXiv:1205.2224v3. Read in full (PDF and TeX).
- [CHM99] M. W. Choptuik, E. W. Hirschmann, R. L. Marsa, "New critical behavior in Einstein-Yang-Mills collapse", PRD 60, 124011 (1999),
  arXiv:gr-qc/9903081. Verified against INSPIRE (record 497094); not read here. It is the issue's reference for the "type III" analogy.

**Update from #223 (range and recount).** The counts above search W ∈ [0.01, 8]. #223 recounted the zero-charge members c = −10, −3, −1.1, −1, −0.9, −0.3, −0.01 with the same code on W ∈ [1e-5, 8] and, for c ≤ −1 (no throat in the static region), with a range-free Sturm count: the same counts (1 mode, 0 at c = −1), W equal to the table to 2.8e-9 (R1) and 3.6e-10 (R2). The other 13 members were not recounted. With a Maxwell field the c = −1 member is not stable: along the locus that keeps the horizon at the area minimum a growing mode with W = 3.8e-3 at q = 0.1 appears (numerical evidence). *Added by the coordinator:* #217's double-null runs of record (`evolution.md`) grow at a rate near κ at c = −1 itself and show a lifetime law with slope 1/κ there, while at c = −0.3 and −1.5 the slope is 1/W; #232 compares the c = −1 case with the master equation.

## Coordinator's review

Under review: `wip/i212-black-universe` at 08a414d, raw `I212-rpi5` at c50bd64. Review folder
`negative-mass-raw/I212-mac-review`, whose README holds the criteria, committed before the first run of `ind212.py`. The
changed files are those the brief names.

- **Rerun on the Mac.**
  - `check`'s verdict lines equal the Pi's.
  - `planted` gives the same verdicts, with W²a² moving by 1.3e-14.
  - `controls`: 14 of 14 fire.
  - The 40 member runs give the same counts, with every W equal to the Pi's to 1.3e-10.
  - `spectral`: 40 of 40. Its first rerun crashed on the review's own command, which omitted the results directory;
    it is kept.
- **ind212**, runs 1 and 2, 3 of 3. It uses a third method: finite differences of u_xx = (V̂ + W²)u on a uniform tortoise
  grid, with Dirichlet ends at ±400 and Richardson extrapolation over two grids. It covers the nine members with no
  throat in the static region.
  - ρ_h and κ equal the table to 1.8e-14.
  - c = −10 … −1.0001 have one negative eigenvalue each, with W equal to the shooting result to 2.4e-10 or better.
  - c = −1 has none.
  - W/κ = 0.99997 at c = −1.0001.
  - It shares the background and V̂ with this note. Those were checked against [BKZ12] (40) and (41) in sympy here,
    and are not re-derived by the review.
  - Mutations: `short-X` fails exactly its target. `wrong-V` fails K2 but is silent on K3, the near-horizon limit.
    `wrong-m` crashes the program instead of failing cleanly (amendment A1 fixed one crash and a second remained). It
    is reported as such.
- *Added by the coordinator:* the throat members (c > −1, eleven of the twenty) are checked here only by this note's two
  methods. The spectral one is weaker near the throat (4.5e-3), so their counts rest mainly on the shooting with the
  complex-path crossing.
- The discrepancies with the printed [BKZ12] (36)–(38) and [BF06] (8) are reported in §2 and §1 and not corrected. They
  belong in `notes/literature.md` with #211.
- `critical_ghost.md` §1.1's row carries an update from this issue.
- **Follow-ups.** #217 (the evolution of a one-mode member, task 2) and #218 (the polar and M-AdS counts, tasks 1 and 3).


# Sourced black bounce and magnetic black universes (#223)

November 2026, #223 (branch `wip/i223-black-bounce`; raw `I223-rpi5` in negative-mass-raw). Code `bounce.py` (it imports the integrator
helpers of `stability.py`), output `bounce_output.txt`. Method and template: #212 (above). Status words as in AGENTS.md: every count
below is *numerical evidence*; the algebra is *checked symbolically* (sympy, inside `bounce.py check`). The result is stated for the two
sources below only: [BW22] say that the stability of an SV geometry may depend on its source.

## Result (#223)

- Of 55 members, 54 have at least one growing radial mode. The exception is B2 at q = 0, which is #212's stable member c = −1; no member
  with a charge or a source term is stable. A Maxwell field does not stabilize the c = −1 member: along the locus that keeps its horizon
  at the area minimum, a growing mode with W = 3.8e-3 (q = 0.1) … 0.168 (q = 0.81) appears at every charge tested.
- A magnetic field adds no radial degree of freedom (*checked symbolically*): in the gauge δφ = 0 the linearized system reduces to one
  second-order equation for Y = δR, whose only new term is `−R A Λ_R` in the potential. The conjecture of the issue holds for both sources.
- Zero charge returns #212's table: seven members (c = −10, −3, −1.1, −1, −0.9, −0.3, −0.01), the same counts, W within 2.8e-9 (R1) and
  3.6e-10 (R2) of the table. The Ellis wormhole through the S code gives W²a² = 1.397052434332 (R2), 1.2e-12 from `w2_output.txt`.
- **Sourced black bounce (S, 16 members, m = M/a from 0 to 10):** every member is unstable. Real growing modes at 14 members;
  a **complex-conjugate pair** (oscillatory growth) at m = 0.2 and 0.3, where there is no real mode; and, for 0.55 ≤ m ≤ 3, a
  **second, exponentially weak real mode** next to the one of rate 0.06 … 0.24 (rates in Sec. 5); at m = 5 and 10 the Sturm count is 1.
- **Magnetic black universes (B, 35 members on three slices):** growing modes at every q > 0 on all three; the charge adds a second, weak
  real mode on B1 and B3 whose rate rises with q (Sec. 5), and one complex pair appears at B1, q = 0.77.
- The counts of the main stage agree between R1 and R2 at every member (rates to 1.3e-7, declared 1e-6). **S9 is NOT MET** at 7 members as
  worded (the range-free Sturm count is 2, shooting in the declared range finds 1; Sec. 6). **S6 is met on both grids at 41 of 55 members** and not
  at 14 (Sec. 8).

## 1. The sources

Static, spherically symmetric, `ds² = A dt² − dx²/A − R² dΩ²`, `R = √(x² + b²)`, units b = 1, a phantom scalar with V(φ) and a
magnetic field F_θφ = q sin θ, `F = 2q²/R⁴`, Lagrangian L(F) ([BW22] (10)–(17), [BBS12] (1), (5)–(7)).

| | metric | scalar and field |
|---|---|---|
| **S**, Simpson–Visser with [BW22]'s source at Q = 0 | `A = 1 − 2m/R`, m = M/a, a = 1 | `φ = arctan x`, `V = 4m cos⁵φ/5`, `L = 12m/(5R⁵)`, i.e. `L ∝ F^(5/4)` ([BW22] (19), (21), (22)) |
| **B**, [BBS12] | `A = R² B`, `B′ = 2(p − x + 2q² arctan x)/(1 + x²)²`, `p = 3m − πq²`, B → 0 at +∞ ([BBS12] (15)–(18)) | `φ = arctan x`, V from (9) = the printed (20), Maxwell field with energy density `q²/R⁴` (7) |

S is a wormhole for m < 0.5 (a > 2m), a one-way wormhole with an extremal horizon at the throat for m = 0.5, a black bounce for m > 0.5;
m = 0 is the Ellis wormhole. B has p > 0 (a de Sitter end); q = 0 is #212's family (3m = −2c/π). Not computed: the RN-SV source of [BW22],
the symmetric p = 0 and the p < 0 members of [BBS12], electric fields, other sources of the same metrics.

Checked (`bounce.py check`, 55 checks, residual 0): the closed forms against the field equations; [BW22] (11)–(13), (16), (17), (21), (22) at
Q = 0 against an independent Einstein tensor in their conventions (signature +−−−, `G^ν_μ = −T^ν_μ`: the printed source is right); [BBS12]
(8)–(12), (15) as printed, and the printed V (20) against (9). No printed equation checked here disagrees with its derivation.

## 2. The equation

Linearization of `R_μν − ε∂_μφ∂_νφ − V g_μν − g_μν S_μ = 0` (ε = −1; `S_t = S_x = Λ + rΛ′/2`, `S_θ = Λ`; Λ = L/2 the NED energy density) and of
`ε□φ = V_φ`, in the gauge δφ = 0, e^(Wt) time dependence, `Y = δR` (φ_c = √2 φ is #212's normalization). Checked symbolically:

- `dF = 0` with `F_θφ = f(t, x) sin θ` forces f constant, `∇_μ(L_F F^(μν)) = 0` holds on the perturbed metric, and `F = 2q²/r⁴` exactly:
  the NED stress tensor depends on the metric only through r, so a radial perturbation has no field of its own;
- the tx equation gives δλ, the scalar equation δν′, the θθ equation the master equation; the tt and xx equations then hold identically,
  in the gauge δφ = 0 and in the gauge δR = 0, and the gauge map `F = −φ₀′Y/R′` takes one onto the other (residual 0).

```math
Y'' = \Big(\frac2x - \frac{A'}A\Big)Y' + \Big(\frac1{R^4} - \frac{A'}{A x R^2} - \frac2{A R^2} - \frac{R\,\Lambda_R}A + \frac{W^2}{A^2}\Big)Y,\qquad
u_{x_*x_*} = (\hat V + W^2)u,\quad u = Y/x,
```

`V̂ = V̂_212 − R A Λ_R`, with `V̂_212` of Sec. 2 above and `R Λ_R = −6m/R⁵` (S), `−4q²/R⁴` (B), so the new term is `+6mA/R⁵` and `+4q²A/R⁴`. Only Λ_R
enters (not Λ, not Λ_RR). The NED extension of [BKZ12] (41), derived here: `A·Ũ = ε φ′²/R′² (R²(V + Λ) − 1) + 2φ′R V_φ/R′ + εV_φφ`, i.e. V → V + Λ in the
first term only; `V̂ = A²Ũ + R_(x*x*)/R` to 1e-29 at four members. `V̂` is real, so growing modes (ω = iW) with real W² are the negative eigenvalues of
`−∂²_(x*) + V̂`; a throat in the static region makes the problem a continuation through x = 0 that is not self-adjoint, and W² may be complex (Sec. 7).

## 3. Method

As #212 (Sec. 3), for the static region outside the outermost horizon x_h, open to x → +∞: ingoing at the horizon, decaying at infinity, Y continued
along a half circle in the complex x-plane through the throat when it lies in the region; a growing mode is a sign change of the Wronskian mismatch on a log grid of W,
refined by Brent. New here: the extremal horizon of the one-way wormhole (`x_* → −∞` like −2/x, a start at x₀ = clip(0.04W, 1e-4, 0.1)); a horizon exactly at the
throat (B2); the main stage W ∈ [1e-3, 8] (80 points R1, 200 R2, rtol 1e-9, 1e-12) and a **low stage** W ∈ [1e-5, 1e-3] (40, 100 points; outer radius
`max(rmax0, 20/W)`, up to 2e6; for B, A and A′ from a series in 1/x at x ≥ 30, which keeps 1e-16 where the closed form loses 1e-11); the **Sturm count** (S9),
the nodes of the zero-energy solution of `u_(x*x*) = V̂u` regular at the horizon for x_h ≥ 0, which has no W grid and so no range; the spectral check S6 as #212's
(Chebyshev, N = 401 and 801, plain grid for box X ≤ 150, `x_* = X sinh(5ξ)/sinh 5` above, box `X = max(X₀, k/W)` from the smallest shooting rate); and a **complex-W shooting** C1
(a 2-D root of the Wronskian mismatch, analytic normalization, from the spectral λ). Planted: P1 Ellis through S, P2 zero charge, P3 a weak bound state W = 3e-5 of a known potential
through the low stage's machinery (8.8e-11 R1, 2.8e-13 R2), P4 a complex bound state W = 0.3 + 0.2i (2.4e-10, 1.4e-13). Declarations, amendments and the development history, with
every control's target and prediction, are in the raw README (A1–A10, D1–D12).

## 4. The counts

W in units 1/b (b = 1, the throat radius, a for S); Sturm = nodes of the zero-energy solution (x_h ≥ 0 only); "unresolved" = the Sturm count exceeds the roots
found down to 1e-5, so the rate lies below 1e-5. R1 and R2 agree at every row (main stage ≤ 1.3e-7, low stage ≤ 3.2e-6 relative). Full table with 10 digits, nodes, κ,
R1 v R2 in `bounce_output.txt`.

### S, sourced black bounce

| m = M/a | a/m | static region | x_h | κ | real modes, W ∈ [1e-3, 8] | real modes, W ∈ [1e-5, 1e-3] | Sturm | complex pair (Re W, ± Im W) |
|---|---|---|---|---|---|---|---|---|
| 0 | ∞ | wormhole (no horizon) | - | - | 1.18197 | - | - | - |
| 0.1 | 10 | wormhole (no horizon) | - | - | 0.55641; 0.608131 | - | - | - |
| 0.2 | 5 | wormhole (no horizon) | - | - | none | - | - | 0.506993 +0.36810i |
| 0.3 | 3.33333 | wormhole (no horizon) | - | - | none | - | - | 0.398398 +0.38720i |
| 0.4 | 2.5 | wormhole (no horizon) | - | - | 0.0107401 | - | - | - |
| 0.45 | 2.22222 | wormhole (no horizon) | - | - | 0.0942891 | - | - | - |
| 0.5 | 2 | one-way, extremal x_h = 0 | 0 | 0 | 0.202427 | - | 1 | - |
| 0.55 | 1.81818 | black bounce | 0.458258 | 0.189363 | 0.00456251; 0.2374 | - | 2 | - |
| 0.6 | 1.66667 | black bounce | 0.663325 | 0.230321 | 0.00239336; 0.241409 | - | 2 | - |
| 0.75 | 1.33333 | black bounce | 1.11803 | 0.248452 | 0.22001 | 0.000553715 | 2 | - |
| 1 | 1 | black bounce | 1.73205 | 0.216506 | 0.177454 | 0.000105581 | 2 | - |
| 1.5 | 0.666667 | black bounce | 2.82843 | 0.157135 | 0.123533 | 1.21486e-05 | 2 | - |
| 2 | 0.5 | black bounce | 3.87298 | 0.121031 | 0.093944 | -; 1 (< 1e-5, unresolved) | 2 | - |
| 3 | 0.333333 | black bounce | 5.91608 | 0.0821678 | 0.063231 | -; 1 (< 1e-5, unresolved) | 2 | - |
| 5 | 0.2 | black bounce | 9.94987 | 0.0497494 | 0.0381208 | - | 1 | - |
| 10 | 0.1 | black bounce | 19.975 | 0.0249687 | 0.0190986 | - | 1 | - |

### B1, fixed mass m = 2/π (c = −3 at q = 0); the horizon crosses the throat at q ≈ 0.7595

| q | x_h | κ | p | real modes, W ∈ [1e-3, 8] | real modes, W ∈ [1e-5, 1e-3] | Sturm | complex pair |
|---|---|---|---|---|---|---|---|
| 0 | 1.11873 | 0.351374 | 1.90986 | 0.282775 | - | 1 | - |
| 0.1 | 1.11063 | 0.351274 | 1.87844 | 0.283023 | 0.000172078 | 2 | - |
| 0.2 | 1.0858 | 0.350865 | 1.7842 | 0.283736 | 0.000723425 | 2 | - |
| 0.3 | 1.04242 | 0.349752 | 1.62712 | 0.00176878; 0.284786 | - | 2 | - |
| 0.4 | 0.976725 | 0.347002 | 1.4072 | 0.00355847; 0.285843 | - | 2 | - |
| 0.5 | 0.881041 | 0.34035 | 1.12446 | 0.00665847; 0.286012 | - | 2 | - |
| 0.6 | 0.73732 | 0.323267 | 0.778886 | 0.0126297; 0.282413 | - | 2 | - |
| 0.7 | 0.480991 | 0.267061 | 0.370479 | 0.0286422; 0.260702 | - | 2 | - |
| 0.74 | 0.264398 | 0.194611 | 0.189523 | 0.0502899; 0.224485 | - | 2 | - |
| 0.76 | -0.0147199 | 0.092972 | 0.0952754 | 0.102543; 0.152211 | - | - | - |
| 0.77 | -1.43678 | 0.11172 | 0.047209 | 0.0395523 | - | - | 0.12126 +0.09153i |
| 0.775 | -2.8048 | 0.152503 | 0.0229402 | 0.0737792 | - | - | - |

### B2, horizon at the throat, m = (2/3π)(1 + q² + π²q²/4); q = 0 is c = −1; κ = p

| q | x_h | κ | p | real modes, W ∈ [1e-3, 8] | real modes, W ∈ [1e-5, 1e-3] | Sturm | complex pair |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.63662 | 0.63662 | none | - | 0 | - |
| 0.1 | 0 | 0.627278 | 0.627278 | 0.00376254 | - | 1 | - |
| 0.2 | 0 | 0.599253 | 0.599253 | 0.0150493 | - | 1 | - |
| 0.3 | 0 | 0.552544 | 0.552544 | 0.033276 | - | 1 | - |
| 0.4 | 0 | 0.487152 | 0.487152 | 0.0571826 | - | 1 | - |
| 0.5 | 0 | 0.403076 | 0.403076 | 0.0850395 | - | 1 | - |
| 0.6 | 0 | 0.300316 | 0.300316 | 0.114664 | - | 1 | - |
| 0.7 | 0 | 0.178873 | 0.178873 | 0.143124 | - | 1 | - |
| 0.75 | 0 | 0.111145 | 0.111145 | 0.155606 | - | 1 | - |
| 0.78 | 0 | 0.0682668 | 0.0682668 | 0.162101 | - | 1 | - |
| 0.8 | 0 | 0.0387468 | 0.0387468 | 0.165884 | - | 1 | - |
| 0.81 | 0 | 0.0237065 | 0.0237065 | 0.167585 | - | 1 | - |

### B3, fixed mass m = 0.2/π (c = −0.3 at q = 0), throat in the static region

| q | x_h | κ | p | real modes, W ∈ [1e-3, 8] | real modes, W ∈ [1e-5, 1e-3] | Sturm | complex pair |
|---|---|---|---|---|---|---|---|
| 0 | -0.926754 | 0.6013 | 0.190986 | 0.97008 | - | - | - |
| 0.03 | -0.93831 | 0.598328 | 0.188158 | 0.00153827; 0.970845 | - | - | - |
| 0.06 | -0.974442 | 0.589143 | 0.179676 | 0.00629203; 0.973147 | - | - | - |
| 0.09 | -1.04014 | 0.572862 | 0.165539 | 0.0146504; 0.977003 | - | - | - |
| 0.12 | -1.14626 | 0.547742 | 0.145747 | 0.0272512; 0.982445 | - | - | - |
| 0.15 | -1.31629 | 0.510551 | 0.1203 | 0.0450033; 0.989519 | - | - | - |
| 0.18 | -1.6088 | 0.454889 | 0.0891983 | 0.069091; 0.998294 | - | - | - |
| 0.2 | -1.95501 | 0.40076 | 0.0653222 | 0.0893425; 1.00514 | - | - | - |
| 0.22 | -2.64605 | 0.320925 | 0.0389328 | 0.113382; 1.01284 | - | - | - |
| 0.23 | -3.39754 | 0.262008 | 0.0247957 | 0.126873; 1.01704 | - | - | - |
| 0.24 | -5.49002 | 0.171477 | 0.0100302 | 0.141271; 1.02149 | - | - | - |

## 5. The second, weak real mode and the charge dependence, as observed

No reading of these modes is offered here.

- **S.** The real root of largest rate: 1.1820 (m = 0, Ellis), 0.5564 and 0.6081 (0.1), none (0.2, 0.3: a complex pair), 0.0107 (0.4), 0.0943 (0.45), 0.2024 (0.5, extremal), 0.2374 (0.55),
  0.2414 (0.6), 0.2200 (0.75), 0.1775 (1), 0.1235 (1.5), 0.0939 (2), 0.0632 (3), 0.0381 (5), 0.0191 (10). For 0.55 ≤ m ≤ 3 there is a second real root, W = 4.56e-3 (m = 0.55, one node), 2.39e-3 (0.6, one node),
  5.54e-4 (0.75), 1.06e-4 (1), 1.21e-5 (1.5); at m = 2 and 3 the Sturm count is 2 and no root is found above 1e-5; at m = 5 and 10 the Sturm count is 1.
- **B1** (m = 2/π): the strong mode stays at 0.2828 … 0.2860 for q ≤ 0.5 and falls to 0.2245 at q = 0.74; the weak mode is W = 1.72e-4 (q = 0.1), 7.23e-4 (0.2), 1.77e-3 (0.3), 3.56e-3 (0.4), 6.66e-3 (0.5),
  0.0126 (0.6), 0.0286 (0.7), 0.0503 (0.74): W/q² = 0.0172, 0.0181, 0.0197, 0.0222, 0.0266. Past the horizon–throat crossing: two modes at q = 0.76 (0.1025, 0.1522), one at 0.77 (0.0396, and the complex pair) and 0.775 (0.0738).
- **B2:** W = 3.763e-3 (q = 0.1), 0.01505 (0.2), 0.03328 (0.3), 0.05718 (0.4), 0.08504 (0.5), 0.1147 (0.6), 0.1431 (0.7), 0.1556 (0.75), 0.1621 (0.78), 0.1659 (0.8), 0.1676 (0.81): W/q² = 0.3763, 0.3762, 0.3697, 0.3574, 0.3402, …
  one mode at every q > 0 and none at q = 0.
- **B3:** the first mode rises slowly from 0.9701 (q = 0) to 1.0215 (q = 0.24); a second mode appears at every q > 0 tested: W = 1.54e-3 (0.03), 6.29e-3 (0.06), 0.01465 (0.09), 0.02725 (0.12), 0.04500 (0.15),
  0.06909 (0.18), 0.08934 (0.2), 0.1134 (0.22), 0.1269 (0.23), 0.1413 (0.24): W/q² = 1.709, 1.748, 1.809, 1.892, 2.000.

## 6. S9: the Sturm count against the shooting count

The shooting counts of record use the declared range W ∈ [1e-3, 8]; the Sturm count has no range. **S9 is NOT MET as worded at 7 of the 33 members where it applies (x_h ≥ 0)**, as found (not loosened; the low stage was added, A2):

| member | shooting, W ∈ [1e-3, 8] | Sturm | low stage, W ∈ [1e-5, 1e-3] | count |
|---|---|---|---|---|
| S m = 0.75 | 1 (0.220010) | 2 | 5.53715e-4 | 2 |
| S m = 1 | 1 (0.177454) | 2 | 1.05581e-4 | 2 |
| S m = 1.5 | 1 (0.123533) | 2 | 1.21486e-5 | 2 |
| S m = 2 | 1 (0.093944) | 2 | none | 2, one rate < 1e-5, unresolved |
| S m = 3 | 1 (0.063231) | 2 | none | 2, one rate < 1e-5, unresolved |
| B1 q = 0.1 | 1 (0.283023) | 2 | 1.72078e-4 | 2 |
| B1 q = 0.2 | 1 (0.283736) | 2 | 7.23425e-4 | 2 |

Each weak root is found at R1 and R2 (agreement 1e-7 … 3.2e-6) and by the Sturm count. At the 26 other members with x_h ≥ 0 the Sturm count equals the shooting count in range.
At the throat members (no Sturm count) the low stage found no root below 1e-3 at any member; a weak mode there would be seen only by shooting and by the spectral check.

## 7. The complex pairs

Eigenvalues λ = W² of the spectral problem (S-B grid, mapped, box from W_ref = 0.03) with Re λ > 1e-3 and Im λ ≠ 0, at 3 of the 55 members; each refined by shooting at complex W (C1, R1 and R2,
declared: R1 v R2 to 1e-6 and λ to 1e-3 of the spectral value; found: 1.3e-10 … 4.9e-10 and 9e-13 … 5.7e-10):

| member | W (Re ± i Im) | λ = W² | real modes of the shooting scan |
|---|---|---|---|
| S m = 0.2 | 0.5069929 ± 0.3680959i | 0.1215472 ± 0.3732441i | none |
| S m = 0.3 | 0.3983976 ± 0.3872015i | 0.0087956 ± 0.3085203i | none |
| B1 q = 0.77 | 0.1212596 ± 0.0915307i | 0.0063260 ± 0.0221979i | 1 (0.0396) |

At S m = 0.1 the shooting finds two real roots (0.5564, 0.6081) that are close to each other; at m = 0.4 one (0.0107). The complex search is local to the spectral eigenvalues (no argument-principle count),
so a pair at another member that the spectral grids miss is not excluded.

## 8. Checks and controls

| check | what | result |
|---|---|---|
| S1 | closed forms and printed sources against the field equations (sympy, residual 0); printed (20) | 55 of 55 checks pass |
| S2, S3 | float64 background against mpmath (incl. the low-stage series); horizons | 3.7e-13 (declared 1e-11); 9e-41 |
| S4 | linearization, generic A, R, φ₀, Λ(R): background, tt/xx after the master equation (two gauges), gauge map, (41) with V → V + Λ, `dF = 0`, NED equation | all 0 |
| S5, S8 | family coefficients against the generic equation (S and B, 4 members); V̂ against the δR = 0 gauge | ≤ 8e-31; ≤ 8e-29 |
| P1, P2 | Ellis through S; zero charge against #212 (7 members) | 3.3e-10 / 1.2e-12; ≤ 2.8e-9 / 3.6e-10 |
| P3, P4 | weak bound state W = 3e-5; complex bound state 0.3 + 0.2i | 8.8e-11 / 2.8e-13; 2.4e-10 / 1.4e-13 |
| R | main stage R1 v R2; low stage R1 v R2 | 55 of 55 (≤ 1.3e-7); 55 of 55 (≤ 3.2e-6) |
| S9 | Sturm against shooting in range (x_h ≥ 0) | **NOT MET at 7 members** (Sec. 6) |
| S6 | spectral, S-A and S-B | both grids at 41 of 55 members; not at 14 |
| C1 | complex shooting at the 3 complex eigenvalues | 3 of 3 |

S6, the 14 members: **3 with complex eigenvalues** (S m = 0.2, 0.3, B1 q = 0.77: the spectral check finds λ that the real-W shooting cannot, Sec. 7); **7 with a weak mode** (S m = 0.6 at S-B, B1 q = 0.3, 0.4 (S-B), B2 q = 0.1, B3 q = 0.03, 0.06, 0.09)
where a missed mode or a rate off by 3e-3 … 0.5 is the discretization of a weakly bound state (the energy of such a state is hypersensitive); **4 at S-A only** (B3 q = 0.2 … 0.24, second mode missed at N = 401). At the three weak
throat members B3 q = 0.03, 0.06, 0.09 and B2 q = 0.1 the second mode is confirmed only by R1 and R2 of the shooting (a scratch run at N = 2401 moved the spectral rate at B3 q = 0.06 from 0.0092 to 0.0069 against the shooting 0.0063; not part of the record).

Controls (`bounce.py controls`; each `--mutate=` must make its target fire): 26 of the 27 declared controls fire, **one is silent: `cplx_conj` on P4** (the conjugate of W in the left start vector; the growing solution dominates, the wrong
admixture is suppressed by e^(−2 Re W r₀), as #212's horizon-exponent control). `cplx_hypot` (the real, non-analytic normalization; P4 finds no root at R1) was **added after the runs of record and is not part of the declared set**;
its target, P4, was computed in run 3 and it fires again in run 4. Those that change a count: `eps` (3 modes instead of 1 at Ellis), `throat_wall` (0), `ned_off` (2 → 1 at B3, q = 0.15), `scan_range` (0
against Sturm 2), `sturm_sign` (Sturm 0 against 2), `spec_sign` (0 against 2), `rmax_low` (0 against 1). `ned_off` and `cplx_ned` show that the NED term matters; it is guarded by S1, S4, S5, S8. The agreement controls
(`hz_start`, `ned_off`, `rmax_low`) compare a mutated R1 with an unmutated R2. A control that was mis-built (`cplx_ned` mutated the spectral reference too) is described in the raw README (A10).

## 9. What this does not show

- Only radial (ℓ = 0) perturbations in the static region outside the outermost horizon; polar and axial modes, the region beyond the horizon, M-AdS and p ≤ 0 members, electric fields and the
  RN-SV source are not computed; the counts are for the members listed, nothing between them.
- A mode count rests on a sign-change search of the declared ranges, a range-free Sturm count at x_h ≥ 0 and a spectral count; complex pairs are searched only near the spectral eigenvalues;
  no argument-principle count and no rigorous bound.
- Nothing is said about whether a bounce or a universe forms from collapse; these are linear modes of static solutions of two particular sources.
- [BW22] note that the same geometry may be stable with another source; the S results hold for the source with V(φ) and L(F) above. The charge slices B1–B3 are one-parameter cuts of a two-parameter family.

## 10. What #212's README, `critical_ghost.md` §1.1 and #211 should say

Proposed wording, for the coordinator (files not touched here).

**#212 README, Sec. 4 or 11, one added paragraph:** *Range and recount. The counts above search W ∈ [0.01, 8]. #223 recounted the zero-charge members c = −10, −3, −1.1, −1, −0.9, −0.3, −0.01 with the same
code on W ∈ [1e-5, 8] and, for c ≤ −1 (no throat in the static region), with a range-free Sturm count: the same counts (1 mode, 0 at c = −1), W equal to the table to 2.8e-9 (R1) and 3.6e-10 (R2).
The other 13 members were not recounted. With a Maxwell field (#223) the c = −1 member is not stable: along the locus that keeps the horizon at the area minimum a growing mode with W = 3.8e-3 at q = 0.1 appears
(numerical evidence).*

**`critical_ghost.md` §1.1, the [BKZ12] row, result column, to add:** *With the sources of [BW22] (Simpson–Visser metric, phantom scalar and magnetic NED field) and [BBS12] (Maxwell field): growing radial modes at 54 of
55 members tested (16 black-bounce members m = M/a = 0 … 10; 35 magnetic black universes on three one-parameter slices; 4 zero-charge members); the exception is the c = −1 member of [BKZ12] at zero charge.
A complex-conjugate pair at S m = 0.2, 0.3 and at one magnetic member; a second, exponentially weak real mode at 0.55 ≤ m ≤ 3 (W = 4.6e-3 … 1.2e-5, below 1e-5 at m = 2 and 3) and at every q > 0 tested on B1 and B3 (#223, numerical evidence, these
sources only).*

**#211 / `baby_universes.md` §4, T2:** *Done as #223. The black bounce with the source of [BW22] and the magnetic black universes of [BBS12] are linearly unstable under radial perturbations at every member tested but one
(c = −1, q = 0); a magnetic field adds no radial degree of freedom. Nothing is shown about non-radial modes, the region beyond the horizon, other sources of the same metric, or formation from collapse.*

## Reproduce (#223)

From this directory, with `~/Source/negative-mass-venv314/bin/python`:

```
bounce.py check; bounce.py planted; bounce.py planted_low; bounce.py planted_complex
bounce.py members; bounce.py count <id> --res=R1|R2 [--stage=low] --out=<id>_<res>[_low].json
bounce.py report DIR --complex=complex.json; bounce.py spectral DIR; bounce.py complex S_m0.2 S_m0.3 B1_q0.77 --out=complex.json; bounce.py controls --base=DIR
```

`bounce_output.txt` is these outputs of record, assembled from negative-mass-raw `I223-rpi5`. The 220 count runs took 6.5 h of CPU on three processes at nice 10 (a main-stage R2 run 1.5 to 7 min,
a low-stage R2 run 15 s to 1.5 min), 117 MB at most; `check` about 1 min; `spectral` and the controls take minutes to tens of minutes.

## References (#223)

- [SV19] A. Simpson, M. Visser, "Black-bounce to traversable wormhole", JCAP 02 (2019) 042, arXiv:1812.07114. Entry and verification in `../literature.md`.
- [BW22] K. A. Bronnikov, R. K. Walia, "Field sources for Simpson–Visser spacetimes", PRD 105, 044039, arXiv:2112.13198. Read: Secs. II–IV in full for the equations used (the sources (10)–(22)); entry in `../literature.md`.
- [BBS12] S. V. Bolokhov, K. A. Bronnikov, M. V. Skvortsova, "Magnetic black universes and wormholes with a phantom scalar", CQG 29, 245006, arXiv:1208.4619. Read: Secs. 1–3 for (1)–(20); entry in `../literature.md`.
- [BF06], [BKZ12]: as above. [BW22] and [BBS12] were read from the arXiv PDFs (v1 of each); (20) of [BBS12] was transcribed from the page image.

## Coordinator's review (#223)

Review of branch `wip/i223-black-bounce` at a44d49d (`bounce.py` 2db4b58) and raw `I223-rpi5` at 476bed0, recorded in
negative-mass-raw `I223-mac-review` (criteria 56439f9, runs c7b0bb8). On the Mac, nice 10.

- **Scope.** The branch adds `bounce.py`, `bounce_output.txt` and the #223 section of this README, as the brief allowed.
- **Rerun.**
  - `check` gives 55 of 55, every verdict line equal to the runner's.
  - `planted`, `planted_complex` and `complex` (S m = 0.2, 0.3, B1 q = 0.77) give every printed W to its 10 digits.
  - `count` at R2 for S m = 1 and B1 q = 0.5 gives the roots to 1e-14.
  - The other 218 count runs were not rerun.
- **Independent check.**
  - The method: a finite-difference eigenvalue problem for −∂²_{x*} + V̂ in the tortoise coordinate, the coordinator's own
    code, sharing only V̂ and A with `bounce.py`.
  - The members: five whose static region has no throat (S m = 0.75, 1, 5 and B1 q = 0, 0.5).
  - The result: their real modes with W ≥ 0.005 come back within 2.5e-4 of the shooting (most within 3e-6), with no other
    eigenvalue with W ≥ 1e-3. Both mutations fire.
  - Not checked independently:
    - the complex pairs, which sit at throat members where the problem is not self-adjoint;
    - the weak modes below 1e-3. A 5 % change of V̂ adds or moves modes near W ≈ 0.005, so those modes are sensitive to
      the potential.
- **Follow-ups.** The runner's four recommended tasks are filed as #236 to #238: tasks 1 and 2 together, task 3, and
  task 4. The c = −1 member is also the subject of #232.

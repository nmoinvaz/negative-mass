# Radial stability of the phantom black universes

November 2026, #212 (branch `wip/i212-black-universe`; raw `I212-rpi5` in negative-mass-raw). Code `stability.py`, output
`stability_output.txt`, graph `growth_rate.svg` (`rate_graph.py`). Literature note: #211. Context: `../critical_ghost.md` §1.1.

Status words as in AGENTS.md. Every count below is *numerical evidence* (a converged computation); the algebra is *checked symbolically*
(sympy, inside `stability.py`). Nothing here is a proof.

## Result

- Nineteen of the twenty black universes tested (c from −10 to −0.01 in units b = 1, c ≠ −1) have **exactly one growing radial mode**; the
  member c = −1, whose horizon sits at the area minimum, has **none**. R1 and R2 agree to 1.2e-8 or better at every member, and an
  independent spectral solver gives the same counts.
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
| c > 0 | M-AdS wormhole; not computed here |

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
sign flipped, the same code gives three growing modes at c = 0 (S, `--mutate=eps`).

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

Near the horizon, with the scaling ρ = ρ_h t at fixed t and ρ_h → 0, V̂ in x tends to a fixed well (sympy, S7, S7c; the members c = −1.0001 and
−0.9999 follow it to 1.5e-3 and 1.3e-3, S7b, S7d):

| side | limit of V̂ | mode | W |
|---|---|---|---|
| c → −1 from below (ρ_h > 0) | `−(A₁²/2)/cosh²(A₁x/2)`, the ℓ = 1 Pöschl–Teller well | `u = 1/cosh(A₁x/2)` | κ = A₁/2 |
| c → −1 from above (ρ_h < 0) | `(A₁²/2)/sinh²(A₁x/2)`, with the throat's wall | `u = 1/sinh(A₁x/2)`, `Y = ρu` regular | κ = A₁/2 |

with A₁ = A′(ρ_h). The well stays at a fixed shape while it moves to x → −∞ like (1/A₁) ln|ρ_h|, so the bound state keeps W = κ. At c = −1 exactly
there is no well at finite x: V̂ → 0 smoothly at the horizon (the wall and the horizon cancel), and the operator has no bound state. The
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

All 19 members with c ≠ −1 listed above. Each has one unstable mode, the ground state of its operator, so each could serve as the critical
solution of a threshold of the kind the issue calls "type III" after [CHM99] (a critical solution with one unstable mode). Nothing is claimed here beyond that: whether
a threshold in phantom collapse approaches any of them is not tested.

## 8. `critical_ghost.md` §1.1, the row on [BKZ12]

The row says: "Phantom with a specific potential (M-AdS wormholes, 'black universes') | Unstable, except black universes whose horizon sits at
the area minimum". [BKZ12]'s abstract: "all configurations under study are unstable under spherically symmetric perturbations, except for a
special class of black universes where the event horizon coincides with the minimum of the area function". **The quote is fair.** Two points
need adding for it to be read correctly:
- the "class" is one member of a one-parameter family (c = −1, 3πm = 2b), not an open set;
- it is a knife edge: members at c = −1 ± 10⁻⁴ have a growing mode of rate κ to 3e-5, not a small one.

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

- Only radial (ℓ = 0) perturbations in the static region ρ > ρ_h. The Kantowski–Sachs region beyond the horizon, polar ℓ ≥ 1 modes (BKZ12
  treat axial modes only at c = −1) and the M-AdS wormholes are not computed.
- Members of this one family and one potential only. Counts are for the 20 values of c listed; between them nothing is computed.
- The count rests on a sign-change search on [0.01, 8] with the node and spectral guards; no argument-principle count, no rigorous bound.

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

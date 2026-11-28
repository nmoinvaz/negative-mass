# A compactly supported datum on a full slice of the c = −1 black universe (#258)

Code [full_slice.py](full_slice.py) (constraints, evolution, checks, mutations), output of record [full_slice_output.txt](full_slice_output.txt), assembled from
negative-mass-raw `I258-mac/runs`; the design (approved by the coordinator) and the criteria, written before the runs of record, are in `I258-mac/design.md` and
`I258-mac/README.md`. Builds on #248 ([horizon_sector.md](horizon_sector.md)), #261 ([horizon_lemma.md](horizon_lemma.md)) and #262 (`cminus1_mode.md` on its branch).
Status words as in AGENTS.md: the slice equations and the first integral Θ are *checked symbolically* only where stated, the rest of the algebra is checked against exact
solutions; every number from the solver is *numerical evidence*. No sentence here says a result is new (#260 reads the literature).

## Answer

**c₂ is free in compactly supported data on a full slice, and the datum evolves to the horizon coefficient the slice predicts.** Datum G (δφ and ∂_Nδφ supported in
|s| ≤ 0.75, amplitude 1e-3, δF = 0 on the slice) solves both constraints on the slice V + U = 0 from the left exterior through S to the right exterior, and the solution is
exactly zero outside the support. It has c₂ = −6.653515e-4. The evolved horizon values are linear in V, and c₂ of their fit equals the slice value to 2.7e-8 relative after
Richardson (bar 1e-6). 17 of 17 checks pass; the 15 mutations fire the declared sets (15 of 15).

- **δM_ADM = 0 is an identity, not a condition.** A first integral Θ of the sum constraint equals −c₁ at S and −δM_ADM/κ at infinity (Sec. 3). For compact δφ, Θ vanishes
  outside the support, so the mass mode is never excited. The conditions that matter are two tail conditions: no constant tail of Y (C1) and no velocity tail (C2).
  The mutation that violates δM_ADM = 0 is a left start with a mass (`left_mass`: δM = 1.00e-4, c₂ changes by 3.1e-4).
- **c₁ = Y(S) = 0 for compact data and c₁ = δM_ADM/κ in general** (Sec. 3). The static sector of #232 and #259 is the mass perturbation.
- **Data supported away from S on both sides have c₂ ≠ 0** (O1: c₂ = −3.18e-4 at amplitude 1e-3): the constraint solution fills the gap. One-sided support away from S gives
  c₂ = 0 exactly (E7). The sentence of README.md and `horizon_sector.md` on data supported away from the horizon needs the qualification of Sec. 7.
- The exact growing mode Y = c₂VR of #262, put through the same code path from its own Y, is reproduced at second order (E9, 6.5e-8 at h = 0.00125). It needs M = 0, which
  the constraints do not fix.

```
 left exterior                         S                        right exterior
 ρ → ∞ ····· zero ····· [ support ] ···┼··· [ support ] ····· zero ····· ρ → ∞      slice T = 0:  V = s, U = −s,  Π = s²
 s = −∞                            s = 0                               s = +∞

                  U ↑   future interior
          V = 0    ╲ │ ╱   U = 0 : H⁺, V > 0   (c₂ fitted here)
                    ╲│╱
        ─────────────S─────────────→ V
                    ╱│╲
        evolved: triangle D⁺[0, 3] (right half, V ≤ 3) and diamond D⁺[−0.5, 0.5] (both future horizons, interior)
```

## 1. The slice equations

Chart of #248, metric `−2e^{−F} dV dU + r² dΩ²`, ε = −1, c = −1. On the slice T = ∂_V − ∂_U = d/ds, N = ∂_V + ∂_U, ∂_t = κ s N. Primes on background functions are d/dΠ.
Perturbation (δr, δF, δφ), gauge δF = δF_N = 0; n = ∂_N δr; free data δφ(s), δφ_N(s) = ∂_N δφ (C^∞ bumps b(x) = e·exp(−1/(1−x²))). C_VV ± C_UU of #261's lemma,
reduced with the evolution equation for δr_VU (`linear.py`'s B0, Bp, Bm; hand derivation):

```
H:  δr'' + P δr' + Q δr = −[E 𝒱′ r δφ + ε r s φ′ δφ'] − [E (1/r − 𝒱 r) δF + s r′ δF']        (' = d/ds, E = e^{−F})
M:  n' + s F′ n = −s r′ δF_N − ε r s φ′ δφ_N
P = sF′ + 2sr′/r,   Q = E𝒱 + E/r² − 2s²r′²/r² + ε s²φ′²
c₁ = δr(0),   c₂ = ∂_V δr(S) = ½(δr'(0) + n(0)),   c₂′ = ∂_U δr(S) = ½(n(0) − δr'(0))      (g = ρR/√2 = 0 at S, so Y = δr there)
δm = δr (½ − e^F s² r′²) − r e^F s r′ (s r′ δF + δr'),      δM_ADM = lim_{s→∞} δm
```

Q1 reproduces exact `L_ξ(member)` (cubic ξ, η; c = −1 and −1.5; δF ≠ 0) to 4.3e-11. The first draft of M had a factor ½ on `δr_TN`; the exact-gauge check showed 3e-2, 3e-11 after the fix.
The homogeneous solution of M is `n_h = e^{−F/2}` (∂_t δr = κ s n_h tends to a constant: a velocity tail). H has two, the constant-tail solution `h_c` (δr → const, δm → 0) and the
mass solution `h_M` (δr ∼ ln ρ, δm → finite). The background is time symmetric on T = 0, so H carries δr, δφ and M carries n, δφ_N.

## 2. Data and conditions

Left: (δr, δr', n) = 0 at s = −1.05, outside the support. Beyond the support the solution is homogeneous; it is required to vanish.

| condition | removes | scalar condition on the data | result |
|---|---|---|---|
| C1 | constant tail of Y (δφ = 0 there) | end vector (δr, δr')(1.05) = 0 | the end vectors of all compact data are parallel to one vector (second singular value 1.7e-14): one condition, not two |
| C2 | velocity tail | n(1.05) = 0 | one condition |

Two correction bumps with solved amplitudes meet them: `δφ += k_m b(s/0.6)`, `δφ_N += k_o [b((s−0.3)/0.3) − b((s+0.3)/0.3)]` (G: k_m = −3.2615e-4, k_o = −3.9597e-4).
The solution is then exactly zero outside the support, which is pure gauge in the strongest sense (in the class δF = δF_N = 0 with δφ = 0 a pure-gauge tail has δr = 0).

Parity P: (V, U) → (−V, −U) is an isometry of the background. Data with δφ odd and δφ_N even need neither condition (end vector 1.7e-17) and carry all of c₂ (G: −6.653515e-4). The P-even part of δφ
and the P-odd part of δφ_N need C1 and C2 and contribute nothing to c₂ (all of δr(0), δr'(0), n(0) ≤ 4e-16). The growing mode has this parity (Y odd, Y_N even).

## 3. The first integral Θ: δM_ADM = 0, c₁ = 0 and c₁ = δM_ADM/κ

With μ = s r, a = μP − μ_s (μ_s = r + 2s²r′), c = ε s² r² φ′ and δF = 0, define `Θ = a δr + μ δr' + c δφ`. Then

```
Θ' = μ · H(δr, δφ)         identically
```

Term by term, Θ' = a'δr + (a + μ_s)δr' + μδr'' + c'δφ + cδφ' and a + μ_s = μP, c = μ ε r s φ′. What is left are two conditions on the background:

- `c' = μ E 𝒱′ r`. With c = ε s² r² φ′ and d/ds = 2s d/dΠ this is `(Π r² φ′)_Π = (ε/2) r² E 𝒱′`, the static scalar equation (*proved* by hand from φ_VU = −(r_Vφ_U + r_Uφ_V)/r − (ε/2)e^{−F}𝒱′).
- `a' = μ Q`: the static Einstein equations. Q6 checks both at six points, worst relative 2.5e-9 (*numerical evidence*).

Consequences:

1. For a solution of H, Θ is constant across the support, and on a homogeneous solution (Q7: relative variation 1.9e-11 from s = 1.05 to 3e4).
2. At S, a(0) = −1 and μ(0) = c(0) = 0, so **Θ = −c₁ = −δr(0) everywhere**.
3. At large s, μ δr' = r A δr_ρ/κ and δm → −r A δr_ρ, while a → 0: **Θ(∞) = −δM_ADM/κ** (the ratio Θ/δm is −1.606 and −1.585 at ρ = 16 for the two
   homogeneous solutions, against −1/κ = −1.571: a 1/ρ convergence, not tested beyond ρ = 16). The constant-tail solution has Θ = 0, the mass solution Θ ≠ 0.

So **c₁ = δM_ADM/κ = (π/2) δM_ADM for every solution of H** at c = −1 (proved by hand in #259 for δφ = δF = 0 beyond a finite ρ, [static_c1.md](static_c1.md) Sec. 2: δm = c_b δr − κΘρ/R with c_b = 1/(2ρ²) + O(ρ⁻³)). A left start with zero constants has Θ = 0: c₁ = 0 and δM_R = 0 for
every compact datum, whatever its shape. That is the rank-1 end vector: the mass direction is never excited, only the constant-tail direction, which C1 removes. The physical reading is
that static solutions are critical points of M_ADM. Checks: Q5 (c₁ = 1.4e-15 and −7.5e-17 from both starts), Q4 (δm ≤ 9e-16 out to ρ = 16), `left_mass` (δM = 1e-4 and 2e-4: c₁ = 1.573913e-4
and 3.147827e-4; κc₁ = 1.001984e-4 against δm(3e4) = 1.000000e-4, the same 0.2% as Θ/δm at ρ = 16 (#259: δm(ρ = 1e8) = κ c₁ to 1.1e-13; the 0.2% is the finite-ρ value of δm); and Θ(1.05) = −c₁ = −1.573913478e-4 to 9 digits, `table`).

The static term M of #232 (u_pm + V̂u/2, u = Y/ρ; zero on the exact mode) follows it. M(s) ≤ 1.1e-14 on G and ≤ 6e-17 on G without C1, and −1.286e-3 at s = ±0.2 for δM = 1e-4 (`table`).
So M ≠ 0 exactly when δM ≠ 0 (hence c₁ ≠ 0): the corner-constant classes A0, A1L of #248 carry a mass perturbation.

## 4. Evolution

`horizon_sector._cellK` (the cell of `linear.py`) on two diagonal layers. Layer 0: values at s = jh. Layer 1: the points (h/2)(1, 1) off the midpoints, `U + (h/2)U_N + (h²/8)U_NN`,
`U_NN = U_TT + 4(B0U + BpU_V + BmU_U)`, with U_T, U_TT from the ODEs and the analytic bumps. Error O(h³) per node, O(h²) overall. Layers ℓ ≥ 2 from ℓ − 2 and ℓ − 1.
Right half: V ≤ 3, bounded by H⁺. Diamond: |s| ≤ 0.5, both future horizons and the interior (Π ≥ −0.25; the interior ends at −0.37). The solution beyond s = 1.05 is set to zero, not extrapolated
(extrapolating the ODE's dense output gave a fake floor of 1e-6 in a development probe). h = 0.01 … 0.00125, 138 MB at the finest, 2 s per march.

## 5. Results

17 of 17 pass (`full_slice_output.txt`).

| id | check | bar | result |
|---|---|---|---|
| Q1 | slice ODEs on exact pure gauge, δF ≠ 0 | ≤ 1e-9 | 1.5e-11, 4.3e-11 (c = −1); 3.5e-12, 6.5e-12 (−1.5) |
| Q2 | δm = ε/2 on Schwarzschild scaling + gauge | ≤ 1e-10 | 5.7e-14, 6.0e-14, 2.6e-11 |
| Q3 | tails of G, 1.05 ≤ \|s\| ≤ 3e4 | ≤ 1e-12 | 2.5e-15 (right), 6.9e-17 (left, right start) |
| Q4 | δM_ADM: \|δm\| at 10 ≤ \|s\| ≤ 3e4 | ≤ 1e-12 | 9.0e-16 |
| Q5 | left start against right start: c₂, c₁ | ≤ 1e-12 | c₂ equal to 10 digits; c₁ 1.4e-15, −7.5e-17 |
| Q6 | Θ′ = μH: a′ = μQ, c′ = μE𝒱′r | ≤ 1e-7 | 2.5e-9 |
| Q7 | Θ conserved on homogeneous solutions; zero past the support | ≤ 1e-9; ≤ 1e-12 | 1.9e-11; 1.4e-15, 5.6e-17 |
| E1 | planted evolution, exact pure gauge, error/scale at h = 0.00125; last ratio | ≤ 3e-5; 4.0 ± 0.3 | 6.5e-6 (−1), 1.4e-6 (−1.5); 4.00, 4.00 |
| E2 | planted c₂: c = −1.5 known −1.12925016e-2 (60-digit quadrature); c = −1 seed | ≤ 1e-5 | slice 5.9e-11, fit 1.0e-8; seed 4.5e-6, ratio 4.00 |
| E3 | horizon values linear in V, h = 0.00125; last ratio | ≤ 1e-5; 4.0 ± 0.3 | 4.8e-7; 4.00 (3.1e-5, 7.8e-6, 1.9e-6, 4.8e-7) |
| E4 | G: c₂ of the slice against the evolved fit, relative, Rich(0.005, 0.0025) | ≤ 1e-6 | 2.7e-8 (2.6e-7 at h = 0.01, 3.7e-9 at the next pair) |
| E5 | constraints of the evolved fields, 8 fixed points, \|C\|/term at h = 0.0025; last ratio | ≤ 1e-3; 4.0 ± 0.3 | 9.5e-6; 4.00, 4.00 (\|C\| 7.1e-7, 1.8e-7, 4.4e-8) |
| E6 | diamond: c₂ and c₂′ fits against the slice | ≤ 1e-6 | 5.5e-8, 1.4e-7 |
| E7 | one-sided R: c₂ slice; fit at h = 0.00125; ratio; Y on \|s\| ≤ 0.5 | 0; ≤ 1e-7; 4.0 ± 0.3; 0 | 0; 2.3e-9; 4.00; 0 |
| O1 | two-sided: c₂ slice; fit error; Y in the gap | ≥ 1e-6; ≤ 1e-6; ≥ 1e-6 | −3.18e-4; 5.0e-7; 1.26e-4 |
| E8 | 20 random compact data: \|Rich − slice\|; min \|c₂\| | ≤ 1e-6; ≥ 1e-5 | 1.1e-9; 3.7e-5 (median 3.9e-4, max 9.4e-4) |
| E9 | planted mode Y = c₂VR on [0.3, 1.6]: error/max at h = 0.00125; ratio; M | ≤ 1e-5; 4.0 ± 0.3; ≤ 1e-9 | 6.5e-8; 4.00; 2.8e-11 |

Datum G, the numbers of the Answer: c₁ = 1.4e-15, c₂ = −6.653514770e-4, c₂′ = 3.568161282e-4, Y(±0.3) = −2.566e-4, 1.920e-4 (the support is |s| ≤ 0.75).
Raw fits of c₂ at h = 0.01 … 0.00125: −6.651343e-4, −6.652971e-4, −6.653379e-4, −6.653481e-4; differences from the slice fall by 4.0 per halving.

**What is shared and what is not.** The slice ODEs are a hand reduction of #261's C_VV, C_UU; they share nothing with the march except the background jets (`horizon_sector.jetsK`), the potential
and the linear system's matrices (also used by the layer-1 formula). E5 evaluates the lemma's C_VV, C_UU by differences on the evolved fields and shares none of the ODE reduction. E9's exact
Y comes from #262's closed form, which shares no code with the march. The coordinator's transport check (`ind248.py`) was not rerun on this datum.

## 6. Mutations

15 of 15 fired sets equal the declared ones (`full_slice.py mutations`). Declared before the runs of record from a development run; the target check is in every set.

| name | change | target | fired set | effect |
|---|---|---|---|---|
| H_P | drop 2sr′/r in P | Q1 | Q1 Q3 Q4 Q5 Q6 Q7 E3 E4 E5 E6 E7 E8 E9 O1 | Q1 0.20 (δr, c = −1); E4 2.3e-2; E5 6.1e-3 and no fall with h |
| H_eps | flip ε in Q | Q1 | same | Q1 7.8e-2; E4 2.9e-2; E5 4.6e-3 |
| M_phi | flip the δφ_N source of M | Q1 | Q1 E3 E4 E5 E6 E7 E8 E9 O1 | Q1 1.5 (n); E4 0.68; E5 3.1e-2 |
| M_F | flip sF′ in M | Q1 | Q1 Q3 E3 E4 E5 E6 E7 E8 E9 O1 | Q1 3.4 (n); E4 8.0e-2; E5 6.2e-3 |
| fn_aT | c₂ = ½ δr'(0) | E2, E4 | E2 E4 E6 E8 O1 | E4 0.30 (c₂ −5.111e-4 for −6.654e-4) |
| fn_sign | c₂ = ½(δr'(0) − n(0)) | E2, E4 | E2 E4 E6 E8 O1 | E4 0.86 (c₂ −3.568e-4) |
| l1_sign | layer 1 with −(h/2)U_N | E1 | E1 E2 E4 E6 E8 E9 O1 | E1 error/scale 1.2, ratio 1.00; E4 0.46 |
| l1_noquad | layer 1 without the h² term | E1 | E1 E2 E3 E4 E5 E6 E7 E8 E9 O1 | E1 1.1e-3, ratio 2.01 (first order); E4 5.8e-4, 2.9e-4, 1.4e-4 |
| dm_noT | δm without the (δF, δr') group | Q2 | Q2 | 1.9e-2, 1.9e-1, 1.2 against ε/2 = 5e-4 |
| mass_off | C1 not imposed | Q3, Q5 | Q3 Q4 Q5 E3 E4 | tails 4.2e-4 at both ends; δm ≤ 1.1e-5; c₂ −4.785e-4 (left start), −8.522e-4 (right start) |
| odd_off | C2 not imposed | Q3, Q5 | Q3 Q5 E3 E4 | κ s n → 8.5e-5; c₂ −6.056e-4 against −7.251e-4 |
| asym_S | zero constants at S, integrate outward | Q3 | Q3 Q4 E3 E4 E6 | tails 5.7e-4; δm ≤ 1.6e-5 (1.1e-6 at 3e4); c₂ = 0 for −6.654e-4 |
| left_mass | start at s = −1.05 with a homogeneous solution of δm(3e4) = 1e-4 | Q4, Q5, Q7 | Q3 Q4 Q5 Q7 E3 E4 | δM_R = 1.000e-4; c₁ = 1.6e-4; Θ(1.05) = 1.6e-4; c₂ −3.590e-4 against −6.654e-4 |
| kg_eps | flip ε in c of Θ | Q6 | Q6 | the scalar equation fails |
| mode_M | E9 with δφ(0.3) = 0 (M ≠ 0) | E9 | E9 | Y off by 3.2e-2 at every h; M reaches 3.6e-4 |

Acceptance 3 of the issue: `asym_S`, `mass_off`, `odd_off` and `left_mass` change c₂ (or leave c₂ ambiguous between the two starts) and `left_mass` violates δM_ADM = 0. `mass_off` leaves a constant Y tail whose
δm decays (1.1e-5 at s = 10, 8.1e-7 at 3e4): it violates "decays to pure gauge", not δM_ADM in the limit. The extra fires have reasons: an ODE mutation breaks the rank-1 end vector (Q3 to Q5); a datum with a tail is cut at s = 1.05 by the
evolution (E3, E4); with c₂ = 0 the relative errors of E4, E6 are degenerate.

## 7. Proposed text (the files are not touched here)

**`horizon_sector.md`**, "What this does not show", first bullet, replace by: *A full Cauchy slice to spatial infinity: built in #258 ([full_slice.md](full_slice.md)). Compactly supported δφ and ∂_Nδφ on the slice V + U = 0 through S
solve both constraints with a solution that is exactly zero outside the support (δM_ADM = 0, no tail of Y), and c₂ from the slice equals c₂ of the evolved horizon fit to 2.7e-8 (numerical evidence).*
"The data are free", class R sentence, replace *"c₂ depends only on data in the past of the crossing … data supported away from the crossing and from S do not excite it"* by: *c₂ of a characteristic datum
depends only on data in the past of the crossing; on a full slice the constraints also tie it to the data on the other side of S. One-sided data supported away from S have c₂ = 0 (#258 E7); data supported on both sides away from S
have c₂ ≠ 0 in general (#258 O1: −3.18e-4 at amplitude 1e-3), because the constraint solution is non-zero between the supports.*

**`README.md`**, Result, first bullet: *which data supported away from the horizon crossing satisfy* → *which one-sided data supported away from the horizon crossing satisfy (on a full slice, data supported on both sides of S away from it
have c₂ ≠ 0 in general, #258)*; Sec. 8, last sentence: *which includes data supported away from the horizon* → *which includes one-sided data supported away from the horizon, not two-sided data (#258)*. After the #248 update add: *Update from #258:
a compactly supported datum on a full slice through S exists with c₂ ≠ 0 (c₂ = −6.65e-4 at amplitude 1e-3); its constraint solution is zero outside the support; c₂ of the slice equals c₂ of the evolved horizon fit to 2.7e-8. First-order δM_ADM of compact
data vanishes identically, and c₁ = δM_ADM/κ (numerical evidence; the first integral Θ with a hand-derived scalar-equation condition).*

**Paper IV, `main.tex`** (read only): replace `\pending{\#258: …}` by

```
\paragraph{A compactly supported datum on a full slice (\#258).} On the slice $V+U=0$ through $S$, with $\delta F=\partial_N\delta F=0$, the two linearized constraints reduce to a second-order and a
first-order equation for $\delta r$ and $\partial_N\delta r$ with the data $\delta\varphi$, $\partial_N\delta\varphi$ free. A first integral $\Theta$ of the sum constraint equals $-c_1$ at $S$ and $-\delta M_{\rm ADM}/\kappa$ at infinity, so the
solution generated by compactly supported data from a zero start on the left has $\delta M_{\rm ADM}=0$ and $c_1=0$, and two scalar conditions (no constant tail of $Y$, no velocity tail) make it vanish outside the support. Datum G of
amplitude $10^{-3}$ has $c_2=-6.653515\times10^{-4}$; the evolved horizon values are linear in $V$ and their fit gives $c_2$ to $2.7\times10^{-8}$ (relative; Richardson of $h=0.005,0.0025$). The exact mode $Y=c_2VR$ of \#262, built from its own $Y$ and
evolved through the same code, is reproduced at second order. Data supported on both sides of $S$ away from it still have $c_2\ne0$, since the constraint solution is non-zero between the supports; one-sided data supported away from $S$ have $c_2=0$.
```

and in the sentence "Data with $c_2=0$, among them all data supported away from the crossing" replace *all data supported away from the crossing* by *one-sided data supported away from the crossing*; delete the final sentence of that paragraph
("a compactly supported datum on a full slice … has not been built") and the open item "A compactly supported datum on a full slice to spatial infinity at $c=-1$ (\#258)". The remaining gaps are in Sec. 8.

## 8. What this does not show

- The late-time growth. The march reaches V = 3 (right half) and |s| ≤ 0.5 (diamond); nothing here follows the datum to V ∼ e^{κt} or to I⁺ of the interior. The rate κ is #232's.
- Non-compact data. The class is strict: exactly zero outside the support. Data with decaying tails (for example the profile of Y_∞ of #262, e^{−κx}) satisfy "decays to pure gauge" in a weaker sense and were not built.
- That G is typical. 20 random compact data have c₂ between 3.7e-5 and 9.4e-4; the space of compact data was not classified.
- A planted exact solution through S. E9 is on [0.3, 1.6]: the equation for δφ has g = ρR/√2 ∼ s² in front of its highest derivative, singular at S.
- The initial values of E9's δφ. The homogeneous solutions with Y = 0 are not all pure gauge: one carries δM. With zero initial values (`mode_M`) M reaches 3.6e-4 and the evolved Y is off by 3.2e-2 at every h; the fix sets δφ(0.3) = α·1e-3, α = 1.171509, so that M = 0.
- C1 and C2 are required for compactness, not shown to be required by the physics; the first integral Θ is checked numerically (a′ = μQ), not derived from the Einstein equations here.
- Nonlinear fate; polar modes; other members.
- *Closed in #259:* Θ(∞) = −δM_ADM/κ is proved for data with δφ = δF = 0 beyond a finite ρ; the ratio at ρ = 16 is the finite-ρ value of δm.

## 9. Recommended tasks

1. *Done in #259.* **c₁ and the mass in #259.** *What:* put the slice data of #248's classes A0, A1L (corner constants) through Θ and read their δM_ADM. *Why:* #259's static sector; c₁ = δM_ADM/κ would make A0, A1L mass perturbations. *Acceptance:* κ c₁ equals
   −κΘ(∞) of the class to 1e-6 (Θ at ρ = 16 and the large-ρ fit); a class with c₁ = 0 gives Θ = 0 to 1e-12. *Depends on:* #259, #258.
2. **The exact mode through S.** *What:* an inverse construction of (δφ, δφ_N) from Y = c₂VR with a Frobenius start at s = 0 (g ∼ s²). *Why:* gives a planted case with c₂ ≠ 0 at the crossing. *Acceptance:* c₂ from the slice equals the mode's to 1e-6; M = 0 to 1e-12;
   a wrong initial value (M ≠ 0) fails. *Depends on:* #262, #258.
3. **Late-time growth of datum G.** *What:* move G's evolved data on two null rays into the (p, σ) chart of #232 and evolve to t = 20. *Why:* the growth at κ, from compact data. *Acceptance:* Y(x = 0, t) rate κ to 1e-3; the class R datum decays. *Depends on:* #258, #232.
4. **Non-compact decaying data.** *What:* data with e^{−κx} tails, C1 and C2 relaxed to decay. *Why:* the weaker reading of "decays to pure gauge". *Acceptance:* c₂ against the evolved fit to 1e-6; δM_ADM = 0 by Θ. *Depends on:* #258.

## 10. Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread, about 40 to 60 s per `check`, 12 minutes for the 15 mutations:

```
full_slice.py check [--only=Q1,..] [--mutate=H_P|H_eps|M_phi|M_F|fn_aT|fn_sign|l1_sign|l1_noquad|dm_noT|mass_off|odd_off|asym_S|left_mass|kg_eps]
full_slice.py table
full_slice.py mutations DIR PREFIX            # fired sets of stored logs DIR/PREFIX<name>.log against the declared ones
```

`full_slice_output.txt` is `check`, `table` and `mutations` of the runs of record (`I258-mac/runs`), in that order. `check.2.log` and `table.2.log` are at code commit 7b7805c; the logs of the 14 older mutations are of run 1
(4f31b3f), which differs from 7b7805c only by the mutation `mode_M` and one `table` line, and `check.log` equals `check.2.log` apart from times. Memory 140 MB at most.

## Coordinator's review

Review of branch `wip/i258-full-slice` at c1276a4 and raw `I258-mac` at a1b39db5, recorded in negative-mass-raw `I258-mac-review`.

- **Design.** The design was approved, with the strict compact reading, the 1e-6 bar relative to c₂, E6 and X20 in, and the two-sided
  finding made a declared check (O1 against E7). It was built as approved.
- **Scope.** The branch adds the script, its output and this note.
- **Rerun.** `full_slice.py check` on the Mac gives 17 of 17, every verdict line as recorded.
- **No separate independent code.** The end-to-end check E4 compares the slice functional with the evolved horizon fit, two paths of this
  note's code. E9 plants #262's exact mode, Y = c₂VR, through the same path to 6.5e-8, and that mode is confirmed independently by the
  coordinator (`I262-mac-review/ind262c.py`). The Θ identity (δM_ADM = 0 for compact data, c₁ = δM_ADM/κ) is numerical evidence, with a
  scalar-equation condition derived by hand.
- **Applied.**
  - `horizon_sector.md` (the class R sentence qualified, the slice gap closed) and the README (the one-sided qualification, an update
    from #258) are updated.
  - Paper IV's #258 placeholder is replaced by the datum, with its numbers checked by `paper_numbers.py` (136 of 136). No placeholder
    remains in Paper IV.
- **Follow-ups.**
  - Task 1 (c₁ of #248's classes through Θ) is added to #259.
  - Task 3 (late-time growth of the compact datum) is #279.
  - Tasks 2 and 4 are not filed: E9 already plants the exact mode, and decaying tails are a weaker reading of a check that passed.

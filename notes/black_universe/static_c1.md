# The static coefficient c₁ of the c = −1 black universe's horizon data (#259)

Code [static_c1.py](static_c1.py), output of record [static_c1_output.txt](static_c1_output.txt), assembled from negative-mass-raw `I259-mac/runs`; the criteria, written before the runs of record,
are in `I259-mac/README.md`. Builds on #248 ([horizon_sector.md](horizon_sector.md)), #258 ([full_slice.md](full_slice.md)), #232 ([linear.md](linear.md)) and #269 ([exact_mode.md](exact_mode.md)).
Status words as in AGENTS.md: the algebra of Secs. 1 and 2 is *proved* by hand and *checked symbolically* where stated, every number from a solver is *numerical evidence*. No sentence here says a result is new
(#260 reads the literature; the first-law reading of Sec. 3 was not searched against it).

## Answer

**c₁ = δr(S) is the coefficient of a static mass mode, and κ c₁ = δM_ADM.** Since r_S = 1, δA_S = 8π c₁ is the first-order change of the area of the bifurcation surface, and the identity is the first law δM = κ δA/(8π)
at S. The classes of #248 with c₁ ≠ 0 (A0, A1R, A1L) are mass perturbations; G, L, R are not. The static perturbation c₁ belongs to is exact: r = √L R, F = A + 1/L − 1 (#269, N6; same potential), with c₁ = δL/2,
c₂ = 0 and δM_ADM = δL/π. It is not a neighbouring member δc: the potential of the member is proportional to c. c₁ is read from one point of any slice through S: c₁ = −Θ_gen(s). 4 of 4 symbolic and 11 of 11 numerical
checks pass; the 18 mutations fire the declared sets.

```
   δr(S) = c₁ ──── δA_S = 8π r_S δr(S) ──── δM_ADM = κ δA_S/(8π) = κ c₁
       │                                          ▲
       │  Θ_gen(s) = −c₁ at every s of the slice  │  δm = c_b δr − κΘρ/R  →  −κΘ   (ρ → ∞)
       └──────────────────────────────────────────┘
```

| acceptance of the issue | result |
|---|---|
| A0, A1L: κ c₁ = δM_ADM read through Θ to 1e-6 | Θ_gen(s) + c₁ ≤ 3.2e-8 (A0), 1.1e-8 (A1L), 7.4e-8 (A1R) relative at s = −0.03 … 1.9; δm(ρ = 1e8) of the cut extensions equals κ c₁ to the same 3.2e-8, 1.1e-8, 7.4e-8 (T6) |
| a class with c₁ = 0 gives Θ = 0 to 1e-12 | compact datum G of #258: 1.4e-15 (T5); evolved L, G, R: 2.0e-13, 3.0e-13, 2.9e-13 after Richardson (T7) |
| c₁ as a functional of the slice data, on #232's solutions to 1e-6 | c₁ = δr(p₀) − e^{bp₀} c₂: 5.05e-9 … 4.36e-8 of \|c₁\| for A0, A1R, A1L at p₀ = 0.5 … 6, \|c₁\| ≤ 9.4e-13 for L, R (T9); from one point, c₁ = −Θ_gen at #232's corner (on the slice): 3.7e-8, 1.3e-8, 5.1e-9 (T11) |
| the static perturbation it equals | the exact mode r = √L R through the same code: c₁ = 1/2, δM = 1/π per unit δL, Θ_gen = −c₁ to 1.1e-11 at 14 values of s, through the evolution to 4.9e-9 in c₁ (T4, T8); dM = κ dA/(8π) along the family to 12 digits (T10) |
| each check has a `--mutate=` with a computed effect | 18 of 18 fired sets equal the declared ones (Sec. 6) |

## 1. The first integral with δF

Chart of #248, slice V + U = 0 (V = s = −U, Π = s²; the member's functions depend on Π). Primes on background functions are d/dΠ, on perturbations d/ds; E = e^{−F}, ε = −1. #258 reduced the sum of the VV and UU
constraints to

```
H:  δr″ + P δr′ + Q δr = −[E𝒱′ r δφ + ε r s φ′ δφ′] − [E (1/r − 𝒱 r) δF + s r′ δF′]
```

and found Θ = a δr + μ δr′ + c δφ, μ = s r, a = μP − μ_s, c = ε s² r² φ′, with Θ′ = μ H(δr, δφ) at δF = 0 (a′ = μQ and c′ = μE𝒱′r are #258's Q6, *numerical evidence*). With δF the same algebra leaves
Θ′ = −μ [E (1/r − 𝒱 r) δF + s r′ δF′], and the background identity

```
(μ s r′)′ = μ E (1/r − 𝒱 r)        (T1: 6.5e-11 at c = −1, 2.8e-10 at c = −1.5; S4, checked symbolically: μ s r′ = ρA/(2κ) and the (θθ) equation 𝒱R² = 1 − (ρA)′)
```

turns it into a total derivative: **Θ_gen = Θ + μ s r′ δF is constant on every solution of H.** At S, μ = 0 and a = −r_S = −1, so Θ_gen(0) = −δr(S) = −c₁, and

**c₁ = −Θ_gen(s) at any point s of the slice** (one point of the data: δr, δr′, δφ, δF).

Θ_gen is gauge invariant. On the exact pure-gauge solutions L_ξ(member) (δF ≠ 0, s = −3 … 10, c = −1 and −1.5) it vanishes to 1.1e-9 and 3.5e-11 of its largest term (finite-difference jets, T3), where Θ without the δF term
is 7.3 |δr| at the worst point.

## 2. The exterior in ρ

Beyond a cut s_c (δφ = δF = 0) and with x = ∫dρ/A, H is

```
δr_xx + p δr_x + q δr = 0,   p = −A′/2 + Aρ/R²,   q = −(A/2)(A″ + ρA′/R²),   δr_x = κ s δr′,   R² = 1 + ρ²
```

(T2: p, q equal the chart's P, Q to 1.4e-13 and 7.0e-12). It has the first integral

```
R (δr_x − ½ A′ δr) = κ Θ           (S1, checked symbolically: its ρ-derivative vanishes with q from (vρ) and R²A″ = 2(A − 1))
```

so δr = √A [C + κΘ ∫dρ/(R A^{3/2})]: C is the constant tail (h_c = √A → 1) and the integral grows like ln ρ (h_M). The Misner–Sharp perturbation there is

```
δm = ½(1 − Aρ²/R²) δr − ρ δr_x = c_b δr − κΘ ρ/R,     c_b = ½(1 − Aρ²/R²) − ½ρA′ = 1/(2ρ²) + O(ρ⁻³)      (S2)
```

and R − ρ = 1/(2ρ) + O(ρ⁻³). So δm → −κΘ with a remainder O(ln ρ/ρ²): **δM_ADM = −κΘ(∞) is proved** for data with δφ = δF = 0 beyond a finite ρ, and with Θ = −c₁,
**κ c₁ = δM_ADM** on every full slice of that kind. Θ has the same value at both ends (it is conserved), so δM_L = δM_R = κ c₁. #258 read the mass at ρ = 16, the finite-ρ value of δm; the limit is read here at ρ = 1e8 from
the ODE in u = ln ρ (the closed form of δr agrees with it to 3.1e-13, Θ is conserved along it to 3.1e-13; T2).

Homogeneous starts of #258 with mass (T5): c₁ = 1.573913e-4 and 3.147827e-4 for the starts normalized to δm(3e4) = 1e-4 and 2e-4, κ c₁ = 1.001984e-4 and 2.003969e-4, equal to δm(ρ = 1e8) to 1.1e-13 and 1.2e-13.
δm(ρ = 1e2) = 1.001841e-4 and δm(ρ = 1e4) = 1.001984e-4 for the first. #258's 0.2% was κ c₁ against δm(3e4), the finite-ρ value of a start normalized to δm(3e4) = 1e-4.

## 3. The static mass mode

The exact static solution of #269 (N6) is r = √L R, φ = √2 arctan ρ, F = A + 1/L − 1, a solution of the same field equations with the same potential. Its mass is M(L) = lim m_L = 2 L^{3/2}/(3π)
(m_L = √L (1 + Lρ²(1 − A))/(2R)), so δM_ADM = δL/π per unit δL: the mpmath derivative at ρ = 1e14 is 0.318309886184, 7.8e-15 from 1/π. In the static gauge δr = (δL/2) R, δφ = 0, so Y = (δL/2) R(ρ),
c₁ = Y(0) = δL/2 and c₂ = 0, and κ c₁ = (2/π)(δL/2) = δL/π = δM_ADM.

The horizon of the perturbed solution moves to ρ_h = (π/4) δL and its surface gravity changes, so the regular form is at fixed Π in the Kruskal chart (#248), with q = ∂_L ρ_L(Π):

```
2κΠ q′ = A′(ρ) q − 1 − κ̂ A(ρ),  κ̂ = κ′/κ = A″(0)/A′(0)² = −π²/8,  q(0) = 1/A′(0)          (the homogeneous solution is the scaling gauge ξ = V∂_V + U∂_U)
δr = R/2 + (ρ/R) q,   δφ = √2 q/R²,   δF = κ̂ − dq/dρ,   c₁ = δr(S) = 1/2,   δM_ADM = 1/π     (per unit δL)
```

through the same code as the classes (a series at |Π| ≤ 0.03, an ODE beyond). Θ_gen, δm of the slice code and the exact ∂_L m at fixed Π (T4):

| s | ρ | Θ_gen + c₁ | Θ without the δF term | δm, slice code | δm, exact ∂_L m |
|---|---|---|---|---|---|
| 0.05 | 0.0062 | 1.665e-16 | −4.999536e-01 | 0.2523826708 | 0.2523826708 |
| 0.3 | 0.1906 | −3.997e-15 | −4.628876e-01 | 0.2966908691 | 0.2966908691 |
| 1 | 1.0000 | −2.220e-15 | 1.068003e-01 | 0.3099079532 | 0.3099079532 |
| 5 | 2.9612 | −1.199e-14 | 2.198572e+00 | 0.3416654330 | 0.3416654330 |
| 1e2 | 7.2516 | 3.020e-14 | 7.108162e+00 | 0.3351927028 | 0.3351927028 |
| 1e4 | 14.1884 | 1.990e-13 | 1.509904e+01 | 0.3278660587 | 0.3278660587 |
| 1e6 | 21.2467 | −5.596e-12 | 2.323409e+01 | 0.3248073903 | 0.3248073903 |

At fixed Π the mass δm(ρ) approaches 1/π only as 0.133/ρ ((δm − 1/π) ρ = 0.136823, 0.133897, 0.133262, 0.133150 at ρ = 1e2, 1e3, 1e4, 1e8) because δφ falls as 1/ρ there; Θ_gen is exactly −1/2 at every s.
The cut extensions (δφ = δF = 0 beyond s_c = ±0.5, ±2, ±20, 1e3, with the jump of δr′ of Sec. 1) read δM_ADM(ρ = 1e8) = κ c₁ to 4.2e-13.

At finite L the first law holds along the whole family. With A_H = 4π L R(ρ_h)², ξ = √L ∂_t (unit at infinity, F_L → 1/L) and κ_ξ = √L A′(ρ_h)/2, dM/dL = κ_ξ (dA_H/dL)/(8π) to 12 digits at L = 0.5, 0.9, 1, 1.5, 2, 4
(T10, mpmath at 40 digits). With κ_ξ without the √L the ratio is 1.414 at L = 0.5 and 0.500 at L = 4 (`law_norm`).

**Not a neighbouring member.** The potential of the member c on its profile is V_c = c V₁ (S3, symbolic), so ∂_c is not a perturbation of the fixed theory. Formally c₁(∂_c) = ∂_c R(ρ_h(c)) = R′(0) ρ_h′ = 0.00e+00 (mpmath)
while ∂_c M_ADM = −0.212206592528 (−2/(3π) = −0.212206590789; the line is evaluated at ρ = 1e14 at 40 digits, which limits it to 8e-9), so κ c₁ = 0 against ∂_c M_ADM: the identity does not hold for it, and need not.
The static solutions of the fixed theory near the c = −1 member that are regular at S and have c₁ ≠ 0 are the exact family above.

## 4. #248's classes through Θ_gen

Kruskal solver, rectangle V ∈ [−0.05, 3], U ∈ [−2, 0.06], data of amplitude 1e-3 (corner at (−0.05, −2)), Richardson of h = 0.0025 and 0.00125. Θ_gen is read at s = −0.03, 0, 0.1, 0.25, 0.5, 1, 1.5, 1.9 (4th-order
difference along the slice for δr′); the cut extensions at s_c = 1.0, 1.9 (right) and −0.03 (left) give δM_ADM(ρ = 1e8).

| class | c₁ (fit) | δr(S) | c₂ | κ c₁ | δM, cut 1.0 | δM, cut 1.9 | δM, cut −0.03 | Θ_gen(0.5) + c₁ | Θ_gen(1.9) + c₁ |
|---|---|---|---|---|---|---|---|---|---|
| A0 | 1.6580466e-03 | 1.6580467e-03 | 1.2270433e-02 | 1.0555453e-03 | 1.0555453e-03 | 1.0555453e-03 | 1.0555453e-03 | −5.36e-11 | −5.24e-11 |
| A1R | 5.2367679e-05 | 5.2367683e-05 | 1.0473537e-03 | 3.3338300e-05 | 3.3338302e-05 | 3.3338302e-05 | 3.3338302e-05 | −3.86e-12 | −3.86e-12 |
| A1L | 2.0947073e-03 | 2.0947073e-03 | 1.0071565e-02 | 1.3335321e-03 | 1.3335321e-03 | 1.3335321e-03 | 1.3335321e-03 | −2.33e-11 | −2.31e-11 |
| G | −4.7908676e-12 | 2.7826511e-17 | −9.6108622e-04 | −3.0499610e-12 | −1.8931993e-13 | −9.0780400e-16 | −2.1614802e-14 | −4.89e-12 | −4.79e-12 |
| L | −9.3074865e-13 | −2.5580639e-17 | −9.6108622e-04 | −5.9253299e-13 | −7.4517792e-16 | −5.1261294e-21 | −2.1659138e-14 | −1.00e-12 | −9.31e-13 |
| R | −3.8601142e-12 | 5.2729887e-17 | 5.3538291e-13 | −2.4574250e-12 | −1.8644313e-13 | 2.1250594e-15 | 4.4660737e-17 | −3.89e-12 | −3.86e-12 |

Relative to |c₁|, Θ_gen(s) + c₁ and δr(S) − c₁ are at most 3.2e-8 (A0), 7.4e-8 (A1R), 1.1e-8 (A1L) over the eight values of s, and the c₁ of the horizon fit equals #248's table to 2.3e-7, 1.8e-8, 1.4e-7 (its seven digits).
For G, L, R the horizon fit gives c₁ = −4.8e-12, −9.3e-13, −3.9e-12 (a truncation seed that falls by 4 per halving, ratios 3.95 … 4.01), while the slice values δr(S) are 2.7826511e-17, −2.5580639e-17, 5.2729887e-17 and
max |Θ_gen(s)| = 3.0e-13, 2.0e-13, 2.9e-13 over the eight values of s. The mass of A0 is δM_ADM = 1.0555453e-3 at amplitude 1e-3, on both sides (the left cut at s = −0.03 gives the same value).

The exact mode as ray data through the same solver (T8, amplitude δL = 2e-3, c₁ = 1e-3, δM = 2e-3/π): error over scale 3.17e-02, 7.96e-03, 1.99e-03, 4.98e-04 at h = 0.01 … 0.00125 (ratios 4.00); relative errors
of c₁ of the fit 4.9e-9, δr(S) 3.5e-9, the functional c₁ = δr − V₀ c₂ at V₀ = 0.5, 1: 2.1e-8, Θ_gen 7.6e-9 and δM of the cuts 7.6e-9; c₂ = −3.1e-11. The Richardson values are c₁ = 9.999999951e-04
(exact 1e-3), δr(S) = 9.999999965e-04, δM of the cuts 6.366197700e-04, 6.366197675e-04, 6.366197698e-04 (exact 6.366197724e-04). The largest pointwise error sits at the far corner V = 3, U = 0.05, where the fields are largest.

## 5. c₁ as a functional of the slice data on #232's solutions

#232's chart starts at p = 0, V = 1, so the horizon values it holds are all at V ≥ 1 and the slice through S is not in it. There are two functionals.

**The horizon law extended to V = 0** (T9): c₁ = δr(p₀, H⁺) − e^{bp₀} c₂, c₂ = ½(∂_T δr + ∂_N δr)/(b e^{bp₀}) from the slice p + σ = const crossing H⁺ at p₀ (#248's functional). Linear solver at c = −1, P = 7,
h = 0.005 and 0.0025 (Richardson):

| class | c₁ (horizon fit) | c₂ (fit) | functional at p₀ = 0.5 | 1 | 2 | 4 | 6 | relative deviation at p₀ = 0.5 … 6 |
|---|---|---|---|---|---|---|---|---|
| A0 | 3.034928390e-04 | 9.5897432e-04 | 3.034928278e-04 | 3.034928278e-04 | 3.034928278e-04 | 3.034928273e-04 | 3.034928257e-04 | 3.67e-08 … 4.36e-08 |
| A1R | −1.570796306e-03 | 1.5707963e-03 | −1.570796327e-03 | −1.570796327e-03 | −1.570796327e-03 | −1.570796327e-03 | −1.570796332e-03 | 1.33e-08 … 1.67e-08 |
| A1L | 1.570796335e-03 | 8.9660477e-04 | 1.570796327e-03 | 1.570796327e-03 | 1.570796327e-03 | 1.570796326e-03 | 1.570796326e-03 | 5.05e-09 … 5.79e-09 |
| L | −9.412418480e-13 | −1.0404810e-04 | 4.384784636e-16 | −2.606729213e-15 | −1.066132136e-15 | 3.288240629e-14 | 1.854916683e-13 | |
| R | 3.719453023e-13 | −1.0690158e-14 | 1.070085170e-26 | 2.959211091e-22 | 2.060751343e-16 | −3.862699999e-14 | −1.697389679e-14 | |

The deviations do not grow with the lever arm e^{bp₀} (1.6 at p₀ = 0.5, 220 at p₀ = 6): they are those of the horizon fit, not of the functional. In the Kruskal solver the same functional at V₀ = 0.5, 1 is
within 5.3e-7 (A0), 1.2e-6 (A1R) and 2.4e-7 (A1L) of c₁ (printed in T6, not a criterion: its error is V₀ times that of c₂, and c₂/c₁ = 20 for A1R).

**One point of the t = 0 slice** (T11): #232's corner (p, σ) = (0, 0) is (V, U) = (1, −1) = (s, −s) at s = 1, ρ = 1, and the corner-constant classes have δφ = δF = 0 there. So c₁ = −Θ_gen(1) = −(a δr + μ δr′),
with δr = a₀, δr′ = ∂_V δr − ∂_U δr = a₁R/b − a₁L tanh(bσ*)/b (b = √2 κ, σ* = 16, ∂_V = ∂_p/(bV), ∂_U = ∂_σ/(−Q′(0)), −Q′(0) = b coth(bσ*)), a(1) = −R A′/(2κ), μ(1) = R = √2:

```
c₁ = (R A′/(2κ))|_{ρ=1} a₀ + (a₁L − a₁R)/κ      (tanh(bσ*) = 1 − 2e^{−2bσ*}, within 1e-12 of 1),     R A′/(2κ) = √2 (1 − π/4) = 0.303492827834 at ρ = 1
```

| class | c₁ from Θ_gen at the corner | c₁ of the horizon fit | relative |
|---|---|---|---|
| A0 (a₀ = 1e-3) | 3.034928278e-04 | 3.034928390e-04 | 3.7e-8 |
| A1R (a₁R = 1e-3) | −1.570796327e-03 | −1.570796306e-03 | 1.3e-8 |
| A1L (a₁L = 1e-3) | 1.570796327e-03 | 1.570796335e-03 | 5.1e-9 |
| L, R | 0 | −9.4e-13, 3.7e-13 | |

The ±(π/2) a of A1R and A1L (κ c₁ = ∓a) is this formula. The mass of the corner classes is δM_ADM = (R A′/2)|₁ a₀ + a₁L − a₁R, with R A′/2 = √2 (2/π − ½) at ρ = 1.

## 6. Checks and mutations

15 checks pass in the run of record (`derive` 4 of 4, `check` 11 of 11; 1 s and 150 s).

| id | what | result |
|---|---|---|
| S1 | first integral R (δr_x − A′δr/2) of the ρ-form, symbolic | residual 0 |
| S2 | c_b = 1/(2ρ²) + O(ρ⁻³), R − ρ = 1/(2ρ) + O(ρ⁻³) | exact |
| S3 | V_c/c has no c | d(V/c)/dc = 0 |
| S4 | μ s r_Π = ρA/(2κ) and (μ s r_Π)′ = μE(1/r − 𝒱r) | residual 0 |
| T1 | the identity of S4 on the chart, c = −1, −1.5 | 6.5e-11, 2.8e-10 |
| T2 | ρ-form against the chart; closed-form tail against the ODE; Θ conserved to ρ = 1e8 | 1.4e-13 (p), 7.0e-12 (q), 1.7e-12 (Θ), 2.4e-13 (δm); 3.1e-13; 3.1e-13 |
| T3 | Θ_gen on exact pure gauge, δF ≠ 0 | 1.1e-9, 3.5e-11 (bar 1e-8) |
| T4 | exact static mode: c₁, Θ_gen, δm, δM, κ c₁, seven cuts | 0.5 exactly; 1.1e-11; 7.2e-13; 4.2e-9 (the 0.133/ρ remainder); 7.8e-15; 4.2e-13 |
| T5 | compact datum G; homogeneous mass starts | 1.4e-15; 1.1e-13, 1.2e-13 |
| T6 | A0, A1R, A1L: Θ_gen, δr(S), δM of the cuts, #248's table | 3.2e-8, 7.4e-8, 1.1e-8; 2.3e-7, 1.8e-8, 1.4e-7 |
| T7 | L, G, R | 2.0e-13, 3.0e-13, 2.9e-13; ratios 4.00 |
| T8 | exact mode through the Kruskal solver | 4.98e-4, ratio 4.00; ≤ 2.1e-8 |
| T9 | #232: c₁ functional | ≤ 4.36e-8 (A0), 1.67e-8, 5.79e-9; \|c₁\| ≤ 9.4e-13 (fit), 1.9e-13 (functional) |
| T10 | first law at finite L | ratios 1.000000000000 at six L |
| T11 | c₁ from Θ_gen at #232's corner | 3.7e-8, 1.3e-8, 5.1e-9 |

18 mutations, each with the checks it must fire declared before the run of record (raw README; two of them, `dm_noT` for T5 and `s_pot` for S4, exposed a check that could not fail and were fixed first): in the run of record the
fired sets equal the declared ones (18 of 18; `static_c1_output.txt`). Two blind spots are kept: `fn_c1_noV` and `fn_c1_sign` do not fire T8, because c₂ = 0 for the exact mode and the V₀ term of the functional is not
tested there (T9 tests it); `th_noF` does not fire T5 (δF = 0 on the compact datum).

**What is shared and what is not.** The exact mode (T4, T8) is an exact solution of the field equations; it shares with the code under test the chart's ρ(Π) (through ρ₁ and the table), the closed forms of A, R and the slice code's
dm and Θ_gen; its reference values (δM = 1/π, the first law at finite L) are mpmath computations from A(ρ) alone. The ρ-form and the reading at ρ = 1e8 use closed forms of A, A′, A″ (series in 1/ρ for ρ ≥ 6) and no chart table; the
closed-form tail and the ODE are two methods for one solution. S1–S4 are sympy and share no code. T6–T8 use #248's solver and its slice functional, T9 and T11 use #232's solver, T5 uses #258's constraint ODEs; the quantities read
from them (Θ_gen, δM from the cut, the corner formula) are computed here.

## 7. Proposed text (the files are not touched here)

**`horizon_sector.md`**, the bullet "What c₁ is", replace the second half ("`c₁ = Y(S)` is the static sector … Data with c₁ ≠ 0 exist (classes A0, A1L below).") by: *`c₁ = Y(S) = δr(S)` is the first-order
change of the area of the bifurcation surface, δA_S = 8π r_S c₁ (r_S = 1 at c = −1), and κ c₁ = δM_ADM (#259, [static_c1.md](static_c1.md)): c₁ is the coefficient of the static mass mode, the exact solution
r = √L R, F = A + 1/L − 1 of #269 (c₁ = δL/2, c₂ = 0, δM_ADM = δL/π). It does not grow. The classes A0, A1R, A1L carry mass; G, L, R do not.* "What this does not show", the bullet "The static source M and c₁ …
is not analysed", replace by: *Done in #259: c₁ is the static mass mode, not a neighbouring member δc (the potential of the member is proportional to c); κ c₁ = δM_ADM to 7.4e-8 for A0, A1R, A1L.* Recommended task 3: done in #259.

**`full_slice.md`**, Sec. 3, after "(conditional on the large-s statement in 3)": replace by *(proved by hand for δφ = δF = 0 beyond a finite ρ, [static_c1.md](static_c1.md) Sec. 2: δm = c_b δr − κΘρ/R with
c_b = 1/(2ρ²) + O(ρ⁻³))*; after "κ c₁ = 1.001984e-4 against δm(3e4) = 1.000000e-4, the same 0.2% as Θ/δm at ρ = 16" add *(#259: δm(ρ = 1e8) = 1.001984e-4 = κ c₁ to 1.1e-13; the 0.2% is the finite-ρ value of δm)*.
Sec. 8, the last bullet ("`c₁ = δM_ADM/κ` rests on Θ(∞) = −δM_ADM/κ, whose numerical ratio at ρ = 16 is 1.02 and 1.01 of the asymptotic value."), replace by *Closed in #259: Θ(∞) = −δM_ADM/κ is proved
for data with δφ = δF = 0 beyond a finite ρ; the ratio at ρ = 16 is the finite-ρ value of δm.* Sec. 9, task 1: done in #259.

**Paper IV, `main.tex`** (read only): in the list of open items delete "; the static sector $c_1$ (\#259)", and after the paragraph "A compactly supported datum on a full slice" add

```
\paragraph{The static coefficient (\#259).} $c_1=Y(S)=\delta r(S)$ is the first-order change of the area of $S$, $\delta A_S=8\pi r_Sc_1$, with $r_S=1$. With $\delta F\ne0$ the first integral of the sum constraint is
$\Theta+\mu s\,r_\Pi\,\delta F$; it is constant on the slice, equal to $-c_1$, and vanishes on pure gauge. Beyond the support, with $x=\int\dd\rho/A$, the constraint has the first integral
$R(\delta r_x-\tfrac12A'\delta r)=\kappa\Theta$ and $\delta m=c_b\,\delta r-\kappa\Theta\rho/R$ with $c_b=1/(2\rho^2)+O(\rho^{-3})$, so $\delta M_{\rm ADM}=-\kappa\Theta=\kappa c_1$ (proved by hand, checked
symbolically). The static perturbation is exact: $r=\sqrt L\,R$, $F=A+L^{-1}-1$ (\#269) has $c_1=\delta L/2$, $c_2=0$ and $\delta M_{\rm ADM}=\delta L/\pi$, and along it $\dd M=\kappa_\xi\,\dd A/8\pi$ holds to $10^{-12}$
for $L=0.5,\dots,4$. The classes A0, A1R and A1L of Table~\ref{tab:c2} have $\Theta_{\rm gen}=-c_1$ along the slice to $7.4\times10^{-8}$ and, with the data cut off beyond $s=1$, $1.9$ and $-0.03$,
$\delta M_{\rm ADM}=\kappa c_1$ to the same accuracy; G, L and R have $\Theta=0$ to $3\times10^{-13}$ (numerical evidence). On the solutions of Sec.~\ref{sec:dn} the corner of the chart lies on the slice $V+U=0$, so
for the corner-constant classes $c_1$ is read at one point, $c_1=\kappa^{-1}\big[\tfrac12RA'\,a_0+a_{1L}-a_{1R}\big]_{\rho=1}$, which is $3.034928\times10^{-4}$ for $a_0=10^{-3}$ and $\pm(\pi/2)\,10^{-3}$ for $a_{1L}$,
$a_{1R}$, against the horizon fits to $3.7\times10^{-8}$, $5.1\times10^{-9}$ and $1.3\times10^{-8}$. It is not a neighbouring member: the potential of the member $c$ is proportional to $c$.
```

## Discrepancies and changes made during the work

- `horizon_sector.md` ("What this does not show") says the nature of c₁ ≠ 0 (a δc-like static mode?) is not analysed: it is the static mass mode r = √L R, not a δc.
- `full_slice.md` Sec. 3 states Θ(∞) = −δM_ADM/κ conditionally and reads the ratio at ρ = 16 (1.02, 1.01 of the asymptotic value), and `left_mass` compares κ c₁ with δm(3e4) (0.2%): both are the finite-ρ value of δm; Sec. 2 proves the limit.
  #258's design note said c₁ ∝ δM_L; Θ gives δM_L = δM_R = κ c₁.
- My own bars moved after development runs: T3 from 1e-9 to 1e-8 (finite-difference jets, d2 gave 1.1e-9) and T8's error/scale from 3e-5 to 1e-3 (d1 gave 5e-4); both are named in the raw README before the run of record. A first draft of
  T6 also required the slice functional of the Kruskal solver at V₀ = 0.5, 1 to 1e-6; A1R missed by 1.2e-6 and the requirement was removed before the run of record (the issue asks for the functional on #232's solutions, T9).
- Run 1 (raw `I259-mac/run1`, code 5d0a28b) was stopped when T11 was found; T11 was declared and developed before run 2 (the run of record). The lines of T1–T10 are equal in the two runs (`I259-mac/runs/INFO.txt`).
- Mutations exposed two checks that could not fail, both fixed before the run of record: T5 divided by δM, which `dm_noT` makes negative, and S4 did not exist when `s_pot` was declared.

## 8. What this does not show

- That the first-law reading is the literature's. δA_S = 8π r_S c₁ with r_S = 1 and κ c₁ = δM_ADM are checked numbers; naming them the first law of black-hole mechanics (Bardeen–Carter–Hawking, Iyer–Wald)
  was not checked against a source *(unverified)*, and #260's notes do not cover it.
- *Done in #301* ([static_c1_members.md](static_c1_members.md)): at c = −1.5, c = −0.3 and on Schwarzschild κ r_S c₁ = δM_ADM, with c₁ = R_h/2 + R′_h/(2κ) per unit δL and δM_ADM = −c/π (= K/2, the static constant at the horizon); c₁ is gauge invariant at R′(ρ_h) ≠ 0, where c₂ is pure gauge.
- A proof for tails. Sec. 2 assumes δφ = δF = 0 beyond a finite ρ. The exact mode has δφ ∼ 1/ρ and δF → const at infinity; Θ_gen = −c₁ and δM = κ c₁ hold for it (T4, T10), which is one planted case, not an argument for such tails.
- That the cut extension is a physical datum. It is a hard cut with the jump condition of Sec. 1; it solves the constraints, and by the conservation of Θ_gen every cut gives the same δM (seven cuts of the exact mode, three of each class agree).
- A slice through S for #232's horizon values. Its chart has none, so T9 extends the horizon law to V = 0 (lever arm e^{bp₀}, 220 at p₀ = 6; the observed deviations do not grow with p₀). T11 is exact only because the corner
  lies on V + U = 0; for other data (the pulses) δφ ≠ 0 away from the corner and there is no single-point formula in #232's chart.
- Other static families of the same theory. r = √L R is the one of #269; whether it is the only regular static perturbation with c₁ ≠ 0 was not searched (the space of regular static data at S is (c₁, δφ(S)) before the condition at infinity).
- Anything nonlinear but the exact family (first law at finite L). Polar modes, other potentials.

## 9. Recommended tasks

1. *Done in #301.* **c₁ and the mass at the other members.** *What:* the slice code of this note (Θ_gen, the ρ-form, the exact family r = √L R, which #269 gives for every member) at c = −1.5, c = −0.3 and Schwarzschild, where r_S ≠ 1 and
   R′(ρ_h) ≠ 0. *Why:* the first-law reading δM = κ r_S δr(S) predicts c₁ for every member; only c = −1 was run. *Acceptance:* κ r_S c₁ = δM_ADM of the exact family to 1e-6 at the three members; a wrong r_S (1) fails at c ≠ −1;
   dM = κ dA/(8π) along the family to 1e-12. *Depends on:* #259, #269.
2. **The decaying partner and the basis (c₁, c₂, c₂′).** *What:* put the time reverse of #262's mode, Y = c₂′ U R(ρ(Π)), through the slice code and read c₂′ = ∂_U δr(S) = ½(n − δr′). *Why:* completes the list of regular
   numbers at S (mass, growing, decaying) as exact modes. *Acceptance:* c₂′ of the slice equals the mode's to 1e-6; a random regular datum's (c₁, c₂, c₂′) is reproduced by the three modes to 1e-6; swapping V and U fails.
   *Depends on:* #259, #262.
3. **A conserved quantity for slices that miss S.** *What:* the analogue of Θ_gen for V + U = const ≠ 0 (not time symmetric), so that c₁ is read from one point of a slice like #232's pulses'. *Why:* T11 reads c₁ at one point only
   for the corner classes; the pulses need the horizon-law extrapolation of T9. *Acceptance:* c₁ from one point of the slice equals the fit to 1e-6 on #232's L and G pulses at five points; dropping the ∂_N terms fails.
   *Depends on:* #259, #232.

## 10. Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread:

```
static_c1.py derive [--mutate=s_q|s_pot]                                   # 4 checks, 1 s
static_c1.py check [--only=T1,..,T11] [--mutate=th_noF|th_nophi|th_asign|dm_noT|kappa_off|jump_off|xform_p|xform_q|mode_khat|mode_dF0|fn_c1_noV|fn_c1_aT|fn_c1_sign|law_norm|cor_noa0|cor_dU]   # 11 checks, 2.5 min
static_c1.py table                                                          # the numbers of this note, 5 min
static_c1.py mutations DIR PREFIX                                           # fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones
```

`static_c1_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I259-mac/runs`, code commit fbc5c8f), in that order. Memory 590 MB at most (#232's solver at h = 0.0025).

## Coordinator's review

Merged at 557416b. Review folder `negative-mass-raw/I259-mac-review` at 9f25c3e8, criteria committed (d3d8630d) before the first run of
`ind259.py`. The changed files are those the brief names.

- **Rerun on the Mac** from 38a08d4: `derive` (4 of 4), `check` (11 of 11) and `table` exit 0, and their logs equal the runner's
  apart from timings.
- **ind259**, 3 of 3, from the closed forms of A, R and #269's family alone. The static member's asymptotic mass is
  M(L) = 2L^{3/2}/(3π) in the areal radius, so dM/dL = 1/π at L = 1. c₁ = ∂_L r(S) = 1/2 and κ = 2/π, so κ c₁ = δM_ADM
  symbolically. Along the family dM/dL = κ_ξ (dA_h/dL)/8π holds to 2e-31 at L = 0.7, 1, 1.3, 2, with κ_ξ normalized at
  infinity. The mutations `no-norm`, `c1-one` and `mass-exp` fire.
- Not rechecked independently: T6 to T11 (the evolved classes and #232's solutions). The two bars the runner moved after
  development (T3, T8) are named in its raw README before the run of record.
- Applied: `horizon_sector.md`, `full_slice.md`, and Paper IV (a paragraph on the static coefficient; the open item removed;
  `paper_numbers.py` 142 strings, mutation `c1theta` added).

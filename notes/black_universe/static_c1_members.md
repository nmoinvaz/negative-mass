# κ r_S c₁ = δM_ADM at the black universes c = −1.5, c = −0.3 and on Schwarzschild (#301)

Code [static_c1_members.py](static_c1_members.py), output of record [static_c1_members_output.txt](static_c1_members_output.txt), assembled from negative-mass-raw `I301-mac/runs`; the criteria, written before the
runs of record, are in `I301-mac/README.md`. Builds on #259 ([static_c1.md](static_c1.md)), #261 ([horizon_lemma.md](horizon_lemma.md)), #269 ([exact_mode.md](exact_mode.md)) and #258 ([full_slice.md](full_slice.md)).
Status words as in AGENTS.md: the algebra of Secs. 1 and 2 is *proved* by hand and *checked symbolically* where stated, every number from a solver is *numerical evidence*. No sentence here says a result is new
(the first-law reading was not searched against the literature; #260 does not cover it).

**Conventions.** G = c = 1. r_S = R(ρ_h) is the areal radius of the bifurcation surface S (the length unit is the throat radius of the Ellis profile, R(0) = 1; Schwarzschild r_S = 2M). κ = A′(ρ_h)/2 is the
surface gravity of the Killing field ∂_t normalized by A → 1 at infinity (Schwarzschild κ = 1/(4M) = 1/(2r_S)). δM_ADM is the Misner–Sharp mass at infinity, m = (r/2)(1 − (∇r)²) at ρ → ∞. c₁ = δr(S) is read per unit
of the family parameter: δL of r = √L R, F = A + 1/L − 1 (the Ellis members), δr_S of r_S → r_S(1 + ε) at r_S = 1 (Schwarzschild).

## Answer

**κ r_S c₁ = δM_ADM at every member tested, with c₁ = R_h/2 + R′_h/(2κ) per unit δL.** For the exact family r = √L R of #269 it holds in closed form: κ r_S c₁ = (κ R_h² + ρ_h)/2 = K/2 = −c/π = dM_ADM/dL,
where K = ½(R²A′ − 2ρA) + ρ is the constant of the static equations (#269) evaluated at the horizon, where A = 0. Through #259's slice functional the identity is read at one point of the slice, because a(0) = −r_S
makes Θ_gen = −r_S c₁ at every s, and beyond a cut δM_ADM = −κΘ. The factor r_S is needed: with r_S set to 1 the identity misses by 5.3e-2 at c = −1.5 and 2.7e-1 at c = −0.3. At R′(ρ_h) ≠ 0 the coefficient
c₂ is pure gauge (#261) but c₁ is not: it is gauge invariant at every member. 3 of 3 symbolic and 8 of 8 numerical checks pass; the 20 mutations fire the declared sets.

```
   δr(S) = c₁ ──── δA_S = 8π r_S δr(S) ──── δM_ADM = κ δA_S/(8π) = κ r_S c₁ = K/2 = −c/π
       │                                          ▲
       │  Θ_gen(s) = −r_S c₁ at every s            │  δm = c_b δr − κΘ R′  →  −κΘ   (ρ → ∞)
       └───────────────────────────────────────────┘
```

| acceptance of the issue | result |
|---|---|
| κ r_S c₁ = δM_ADM to 1e-6 at the three members | exact mode through the slice code: c₁(S) against d_L r_S 4.1e-16 (c = −1.5), 1.4e-15 (c = −0.3), 0 (Schwarzschild); Θ_gen + r_S c₁ ≤ 4.5e-12, 8.2e-13, 2.7e-15 at 8–14 values of s; κ r_S c₁ against δM_ADM (mpmath) 1.2e-16, 1.5e-16, 0; the 6–7 cut extensions read at ρ = 1e8 equal δM_ADM to 2.0e-13, 2.7e-13, 1.7e-13 (T4). Through the solver: 9.5e-8 (Θ_gen, c = −0.3) … 8.5e-7 (κ r_S c₁ from the horizon fit against the cuts, c = −0.3), and ≤ 1.9e-7 for A0, A1R, A1L of #248 (T7, T8) |
| a wrong r_S (set to 1) fails at c ≠ −1 | κ·1·c₁ against δM_ADM: 5.31e-2 (c = −1.5), 2.67e-1 (c = −0.3); 0 at c = −1 and at Schwarzschild 2M = 1, where r_S = 1 and nothing can fail; Schwarzschild 2M = 2 (closed forms) 5.0e-1 (T5) |
| dM = κ dA/(8π) along each family to 1e-12 | ratio 1.0000000000000 at L = 0.5, 0.9, 1, 1.5, 2, 4 for c = −1.5, −0.3, −1 and at r_S = 0.5, 1, 2 for Schwarzschild; M_MS = K L^{3/2}/3 to 0 (120 digits) (T2) |
| each check has a `--mutate=` with a computed effect | 20 of 20 fired sets equal the declared ones (Sec. 6) |

## 1. The exact family at a member

Members of #212: R² = 1 + ρ², φ = √2 arctan ρ, A = 1 + cR² + K(ρ + R² arctan ρ), K = −2c/π; the horizon is the simple zero ρ_h of A (A is strictly increasing for c < 0). The static mass mode is #269's N6, r = √L R,
F_L = A + 1/L − 1, φ unchanged, with the same potential. Its horizon moves to A(ρ_h(L)) = 1 − 1/L, so δρ_h = δL/A′(ρ_h) = δL/(2κ), and

```
c₁ = ∂_L [√L R(ρ_h(L))] = R_h/2 + R′_h/(2κ)         (the second term is the horizon moving through the profile; it is zero at c = −1)
κ r_S c₁ = κ R_h²/2 + R_h R′_h/2 = (κ R_h² + ρ_h)/2 = K/2                    (R_h R′_h = ρ_h; K at A = 0)
M_ADM(L) = lim m_L = K L^{3/2}/3 = −2c L^{3/2}/(3π)   so   dM_ADM/dL = K/2 = −c/π
```

(S3, *checked symbolically*: K is constant on R²A″ = 2(A − 1); the identity above; lim ρ(1 − A)/2 = K/3; the mass is read from m_L = √L(1 + Lρ²(1 − A))/(2R)). The mass is the g_tt route as well (by hand, not run): L F_L = 1 − (2K/3)L^{3/2}/r.
The values (T1, mpmath at 40 and 120 digits; the double-precision constants agree to 1e-12):

| member | ρ_h | κ | r_S | R′(ρ_h) | c₁ per δL | δM_ADM per δL | κ r_S c₁ | κ c₁ (r_S = 1) |
|---|---|---|---|---|---|---|---|---|
| c = −1 (#259) | 0 | 0.636619772 | 1 | 0 | 0.5 | 1/π = 0.3183098862 | 0.3183098862 | 0.3183098862 |
| c = −1.5 | 0.339617129 | 0.551681679 | 1.056096489 | 0.321577746 | 0.8195005061 | 3/(2π) = 0.4774648293 | 0.4774648293 | 0.4521034150 |
| c = −0.3 | −0.926753836 | 0.601299801 | 1.363404809 | −0.679734903 | 0.1164811100 | 3/(10π) = 0.0954929659 | 0.0954929659 | 0.0700400683 |
| Schwarzschild, 2M = 1 | 1 (r) | 1/2 | 1 | 1 | 1 (per δr_S) | 1/2 | 1/2 | 1/2 |

At c = −0.3 the two terms of c₁ are 0.6817 and −0.5653: the horizon moves outward in ρ, down the throat's slope (R′ < 0 there), and takes back 83% of the rescaling.
Schwarzschild: m = r_S/2 identically and κ = 1/(2r_S), so κ r_S c₁ = ½ δr_S = δM_ADM for every r_S; K does not exist there.

**Finite L.** With ξ = √L ∂_t unit at infinity, κ_ξ = √L A′(ρ_h(L))/2 and A_H = 4π L R(ρ_h(L))², dM = κ_ξ dA_H/(8π) holds to 1e-13 at six values of L for each member (T2: M is the Misner–Sharp
limit at ρ = 1e20, central differences at 120 digits). Without the √L the ratio is 1.414 at L = 0.5 (`law_norm`); with M ∝ L² it fails (`law_mass`).

## 2. The slice functional at any member, and the gauge of c₁

**Θ_gen.** #259's derivation uses nothing of c = −1 but the background identity (μ s r_Π)′ = μE(1/r − 𝒱r). With r = R(ρ), E = A/(2κ²Π), ρ_Π = A/(2κΠ) and the (θθ) equation 𝒱R² = 1 − ½(A(R²)′)′, it holds
for generic R(ρ) and A(ρ) (S2, *checked symbolically*; this covers Ellis members and Schwarzschild at once). The coefficient of δr is a = r(s²F_Π − 1) = −R A′/(2κ), so **a(0) = −R(ρ_h) = −r_S**, because A′(ρ_h) = 2κ
(S2; numerically a(s = 1e-6) + r_S ≤ 9.3e-12, T3). Hence Θ_gen(0) = −r_S δr(S) and **c₁ = −Θ_gen(s)/r_S at any point s of the slice**. At c = −1 r_S = 1 and this is #259's statement.

**The exterior.** Beyond a cut (δφ = δF = 0) with x = ∫dρ/A, H is y_xx + p y_x + q y = 0 for y = δr, with

```
p = −A′/2 + A R′/R,   q = −(A/2)(A″ + A′R′/R)          (T3: equal to the chart's P, Q to 1.7e-13 and 6.4e-12 at c = −1.5, 1.9e-14 and 2.7e-12 at c = −0.3, 1.7e-15 and 2.6e-14 for Schwarzschild)
R (y_x − ½A′y) = κΘ,      δm = ½(1 − A R′²) y − R R′ y_x = c_b y − κΘ R′,      c_b = ½(1 − A R′² − R R′ A′)
```

(S1, *checked symbolically*: the first integral from the chart's Q, for the Ellis members with generic A and for Schwarzschild.) For the Ellis members A − 1 = −2K/(3ρ) + O(ρ⁻³), R′ = 1 − 1/(2ρ²) + …, so c_b = 1/(2ρ²) + O(ρ⁻³)
(by hand; not in the symbolic checks), y = O(ln ρ), and δm → −κΘ: **δM_ADM = −κΘ(∞)**. For Schwarzschild c_b = 0 exactly and δm = −κΘ at every ρ. With Θ = −r_S c₁: **κ r_S c₁ = δM_ADM** on every full slice whose data vanish beyond a
finite ρ. The slice at c = −0.3 passes the throat R′ = 0 (ρ = 0, s ≈ 0.6) outside the horizon; no coefficient is singular there (T3 compares p, q at s = 0.5 and 1).

**Gauge.** The lemma of #261 (item 2) gives L_ξ r = −V ξᵁ R′(ρ_h) ρ′(0) on U = 0, which vanishes at V = 0 for every ξ regular at S. So δr(S) = c₁ is invariant under every linear gauge transformation at every member; c₂ shifts by
−η(0) R′(ρ_h) ρ′(0) for ξ = η(U)∂_U. T6 runs exact pure gauge L_ξ(member) (η(0) = 0.03, δF ≠ 0): Θ_gen = 0 over its largest term to 1.1e-9 (c = −1), 3.5e-11, 2.4e-10, 6.0e-12 (Schwarzschild); δr(S) = 0 and the
horizon-fit c₁ ≤ 6.5e-17; and c₂ equals −η(0)R′(ρ_h)G(0) (mpmath, G(0) = dρ/dΠ at Π = 0) to 2.3e-13 and 9.3e-14:

| member | R′(ρ_h) | c₂ of the gauge solution (η(0) = 0.03) | c₁ of the gauge solution |
|---|---|---|---|
| c = −1 | 0 | 0 | 0 |
| c = −1.5 | 0.321578 | −1.129250e-2 (#248's N6: −1.12925016e-2) | 4.6e-18 |
| c = −0.3 | −0.679735 | +1.359160e-1 | −6.5e-17 |
| Schwarzschild | 1 | −1.103638e-2 | 2.3e-18 |

So at c ≠ −1 and on Schwarzschild the pair (c₁, c₂) of the horizon law has one invariant, c₁, and one gauge number, c₂; the area of S, δA_S = 8π r_S c₁, is invariant. This does not say what the invariant growth of the
event-horizon area is when c₂ is pure gauge (Sec. 8).

## 3. The exact mode through the slice code (T4)

The fields of the mode at fixed Π (Sec. 3 of #259's note, generalized: q = ∂_Lρ_L(Π) solves A q′ = A′q − 1 − κ̂A, κ̂ = A″(ρ_h)/A′(ρ_h)² = −1/(2κ²R_h²), regular at ρ_h with q(ρ_h) = 1/A′(ρ_h); a Taylor series about ρ_h
from R²A″ = 2(A − 1) for |ρ − ρ_h| ≤ 0.2, an ODE beyond) are put through Θ_gen, dm and the cut extensions. Per unit δL (δr_S):

| member | c₁(S) | Θ_gen + r_S c₁ (max over s) | c₁ = −Θ_gen/r_S | δm slice against exact | δM_ADM exact | cut extensions at ρ = 1e8 |
|---|---|---|---|---|---|---|
| c = −1 | 0.500000000000 | 6.5e-12 | 6.5e-12 | 6.7e-13 | 0.318309886184 | 2.1e-13 |
| c = −1.5 | 0.819500506148 | 4.5e-12 | 4.5e-12 | 7.9e-13 | 0.477464829276 | 2.0e-13 |
| c = −0.3 | 0.116481110020 | 8.2e-13 | 8.2e-13 | 2.9e-13 | 0.095492965855 | 2.7e-13 |
| Schwarzschild | 1.000000000000 | 2.7e-15 | 2.7e-15 | 1.0e-15 | 0.500000000000 | 1.7e-13 |

(relative; s from 0.05 to 1e5 for the Ellis members and 0.05 to 10 for Schwarzschild, whose chart loses Newton convergence beyond Π ≈ 150; cuts at s = ±0.5, ±2, ±20, 1e3 and ±0.5, ±2, ±10.) The tables of `static_c1_members_output.txt` give
Θ_gen + r_S c₁ at each s: for c = −1.5 and −0.3 it is below 4e-12 and 1.4e-13, while Θ without the δF term reaches 32 and 7.7. At fixed Π the mass δm(ρ) approaches δM_ADM like 1/ρ with
(δm − δM)ρ = 0.133150 (c = −1), 0.013529 (c = −1.5), 0.378030 (c = −0.3) at ρ = 1e8, as in #259 Sec. 3: a mode with δφ ∼ 1/ρ is read at a cut, where δφ = 0 beyond.

**What is shared.** The reference values (c₁, δM_ADM, the first law) are mpmath computations from A(ρ) alone at 120 digits (derivatives of √L R(ρ_h(L)) and of lim m_L at ρ = 1e20); the mode's q is an ODE or series in ρ;
the slice code is #258/#259's (Θ_gen, dm) on the chart's coefficients. The ρ-form and its tail use closed forms of A, A′, A″ (series in 1/ρ for ρ ≥ 6) and no chart table; the ODE and the closed-form tail
(y = √A[y₀/√A₀ + κΘ ∫dρ/(R A^{3/2})]) agree to 4.9e-13 (T3).

## 4. Through the solver (T7, T8)

Rectangle V ∈ [−0.05, 3], U ∈ [−2, 0.06], h = 0.005 … 0.000625, Richardson of the two finest; every number is relative to the exact value or to max(|c₁|, 1e-5).

**The exact mode as ray data** (amplitude 2e-3; T7). The pointwise error over scale falls by 4.00 at each halving at every member: 7.96e-3 … 1.25e-4 (c = −1, 4.98e-4 at h = 0.00125 as in #259), 9.90e-5 … 1.55e-6 (c = −1.5),
1.36 … 2.13e-2 (c = −0.3), 6.16e-6 … 9.62e-8 (Schwarzschild).

| member | c₁ of the fit | δr(S) | Θ_gen(s) + r_S c₁ | δM of the cuts against exact | κ r_S c₁(fit) against δM of the cuts | c₂ |
|---|---|---|---|---|---|---|
| c = −1 | 2.6e-10 | 2.2e-10 | 4.3e-10 | 4.3e-10 | 1.8e-10 | −2.0e-12 |
| c = −1.5 | 2.9e-11 | 7.0e-13 | 8.8e-11 | 8.8e-11 | 8.3e-11 | −3.9e-14 |
| c = −0.3 | 9.0e-7 | 5.2e-8 | 9.5e-8 | 9.5e-8 | 8.5e-7 | 1.3e-10 |
| Schwarzschild | 2.3e-12 | 4.9e-14 | 1.3e-10 | 1.3e-10 | 1.2e-10 | −1.6e-15 |

At c = −0.3 the chart has e^{−F} = 10.9 at s = 0.05 (S sits at ρ = −0.93 while the chart is scaled at ρ = 1), the second-order constant is larger, and the horizon-fit c₁ reaches 9.0e-7: 1.1 times
inside the issue's bar. The slice readings (δr(S), Θ_gen, cuts) are 5e-8 … 1e-7 there.

**#248's classes** (amplitude 1e-3, corner (−0.05, −2); T8). Per member, c₁ of the horizon fit, κ r_S c₁, and the largest relative deviation over s = −0.03 … 1.9 of Θ_gen + r_S c₁, over the three cuts of δM − κ r_S c₁, and of δr(S) − c₁:

| member | class | c₁ | κ r_S c₁ | Θ_gen | δM of the cuts | δr(S) | c₂ (fit) |
|---|---|---|---|---|---|---|---|
| c = −1 | A0 | 1.6580467e-3 | 1.0555453e-3 | 2.9e-9 | 2.9e-9 | 2.9e-9 | 1.227e-2 |
| c = −1 | A1R | 5.2367683e-5 | 3.3338302e-5 | 6.6e-9 | 6.6e-9 | 6.4e-9 | 1.047e-3 |
| c = −1 | A1L | 2.0947073e-3 | 1.3335321e-3 | 9.2e-10 | 9.2e-10 | 9.1e-10 | 1.007e-2 |
| c = −1.5 | A0 | 1.1955147e-3 | 6.9654161e-4 | 8.1e-11 | 8.2e-11 | 7.6e-11 | 3.761e-3 |
| c = −1.5 | A1R | 4.8374498e-5 | 2.8184389e-5 | 7.8e-10 | 7.8e-10 | 1.7e-10 | 9.675e-4 |
| c = −1.5 | A1L | 1.9349799e-3 | 1.1273756e-3 | 3.6e-11 | 3.5e-11 | 2.7e-11 | 2.694e-3 |
| c = −0.3 | A0 | 2.9850145e-3 | 2.4471598e-3 | 1.4e-7 | 1.4e-7 | 1.4e-7 | 3.884e-2 |
| c = −0.3 | A1R | 7.8267371e-5 | 6.4164768e-5 | 1.9e-7 | 1.9e-7 | 1.9e-7 | 1.565e-3 |
| c = −0.3 | A1L | 3.1306952e-3 | 2.5665910e-3 | 9.3e-8 | 9.3e-8 | 9.3e-8 | 3.445e-2 |
| Schwarzschild | A0 | 1.0397402e-3 | 5.1987008e-4 | 7.2e-11 | 7.2e-11 | 1.4e-12 | 7.9e-4 |
| Schwarzschild | A1R | 4.8088938e-5 | 2.4044469e-5 | 8.7e-10 | 8.7e-10 | 2.5e-10 | 9.62e-4 |
| Schwarzschild | A1L | 1.9235575e-3 | 9.6177876e-4 | 5.3e-12 | 2.8e-12 | 2.1e-13 | 1e-15 |

The c = −1 rows equal #259's table (1.658047e-3, 5.236768e-5, 2.094707e-3; κ c₁ = 1.0555453e-3, 3.3338302e-5, 1.3335321e-3). At c = −1.5 and −0.3 κ c₁ (r_S = 1) is 6.595e-4 and 1.795e-3 for A0: the factor r_S is what
the cuts see. The pulse class G has c₁ = −5.2e-13, −7.9e-14, −2.5e-12 at c = −1, −1.5, −0.3 (a truncation seed: it falls by 4 per halving at c = −1 in #259), |Θ_gen + r_S c₁| ≤ 3.4e-12 and cuts |δM| ≤ 1.3e-14: a
class with c₁ = 0 gives Θ_gen = 0 at every member. c₂ at Schwarzschild A1L is 1e-15 and 9.6e-4 for A1R: the data of the class fix a gauge, and c₂ is not a number of the solution there.

## 5. Checks and mutations

All 11 pass in the run of record (`derive` 3 of 3, 0 s; `check` 8 of 8, 128 s; `table`).

| id | what | result |
|---|---|---|
| S1 | chart Q → ρ-form q; first integral R(Ay′ − A′y/2), Ellis and Schwarzschild | residual 0 |
| S2 | (μ s r_Π)′ = μE(1/r − 𝒱r), a = −R A′/(2κ), generic R, A | residual 0 |
| S3 | K constant; κ r_S c₁ = K/2; ρ(1 − A)/2 → K/3; Schwarzschild | residual 0 |
| T1 | constants against mpmath; c₁ closed form; δM = −c/π; κ r_S c₁ = δM from the constants | ≤ 1e-12; 5.8e-16 |
| T2 | first law at finite L | ratios 1.0000000000000 |
| T3 | ρ-form, a(0) = −r_S, κΘ and δm in ρ against the s-form, closed-form tail, conservation to ρ = 1e8 | ≤ 1.4e-11 |
| T4 | exact mode through the slice code | ≤ 6.5e-12; 2.7e-13 |
| T5 | wrong r_S | fails by 5.3e-2, 2.7e-1, 5.0e-1; holds to 0 where r_S = 1 |
| T6 | pure gauge | 1.1e-9; c₂ to 2.3e-13 |
| T7 | exact mode through the solver | ratios 4.00; ≤ 9.0e-7 |
| T8 | classes through the solver | ≤ 2.5e-7 (G at c = −0.3, in units of 1e-5); ≤ 1.9e-7 for A0, A1R, A1L |

Twenty mutations, each with the checks it must fire declared before the runs of record (raw README, effects computed in development run d1): in the run of record the fired sets equal the declared ones (20 of 20;
`static_c1_members_output.txt`). `rs_one` (r_S := 1) fires T1, T4, T7, T8; `ctl_true` fires T5 (the control must be able to fail); `c1_noshift` (c₁ = R_h/2) fires T1 with c₁ = 0.528 for 0.8195 at c = −1.5; `th_noF`,
`th_nophi`, `th_asign`, `dm_noT`, `jump_off`, `xform_p`, `xform_q`, `rho_ellis` (Schwarzschild with the Ellis R), `mode_khat` (ratio 1.00: no convergence), `mode_dF0`, `kappa_off`, `law_norm`, `law_mass`, `gauge_noRp`
and the three symbolic ones fire what the README lists. Blind spots kept: `dm_noT` does not fire T3 (the s-form and ρ-form δm share the mutation; T4 sees it), `mode_*` do not fire T8, and `gauge_noRp` is silent where
R′_h = 0.

## 6. Proposed text (the files are not touched here)

**`static_c1.md`**, Sec. 8, replace the bullet "Any member but c = −1. …" by: *Done in #301 ([static_c1_members.md](static_c1_members.md)): at c = −1.5, c = −0.3 and on Schwarzschild κ r_S c₁ = δM_ADM, with
c₁ = R_h/2 + R′_h/(2κ) per unit δL and δM_ADM = −c/π; c₁ is gauge invariant at R′(ρ_h) ≠ 0 while c₂ is pure gauge.* Sec. 9, task 1: done in #301. The Answer's "κ c₁ = δM_ADM" gains "(κ r_S c₁ in general; r_S = 1 here)".

**`horizon_lemma.md`**, in "The members", add the row `| c = −0.3 black universe | −0.679735 | 0.601300 | c₂ pure gauge |` and after the Lemma's item 3 the item: *5. c₁ is invariant. By item 2, L_ξ r = −V ξᵁ R′(ρ_h)ρ′(0) on U = 0
vanishes at V = 0, so c₁ = δr(S) is the same in every gauge at every member, also where R′(ρ_h) ≠ 0 and c₂ is pure gauge. δA_S = 8π r_S c₁ with r_S = R(ρ_h) (#301: Θ_gen = 0 and δr(S) = 0 on exact pure gauge at c = −1, −1.5, −0.3 and on
Schwarzschild; c₂ shifts by −1.12925e-2 at c = −1.5 and +1.35916e-1 at c = −0.3 for η(0) = 0.03).*

**Paper IV, `main.tex`** (read only), in "The static coefficient": replace "$\delta A_S=8\pi r_Sc_1$ with $r_S=1$." by "$\delta A_S=8\pi r_Sc_1$ with $r_S=R(\rho_h)$ ($r_S=1$ at $c=-1$)." and after the sentence ending "so $\kappa c_1=\delta M_{\rm ADM}$,"
insert, at the end of the paragraph (before "So data with $c_1\ne0$ change the mass."):

```
At every member the same family has $c_1=\delta L\,[R_h/2+R'(\rho_h)/(2\kappa)]$, $M_{\rm ADM}=-2cL^{3/2}/(3\pi)$ and $\kappa r_Sc_1=\delta M_{\rm ADM}$ (the constant $K=\tfrac12(R^2A'-2\rho A)+\rho$ of the static
equations at the horizon gives $\kappa r_S^2+\rho_h=K$; checked symbolically), and $\Theta_{\rm gen}(0)=-r_Sc_1$. Through the slice code the identity holds to $4\times10^{-13}$ at $c=-3/2$ ($r_S=1.056$),
$c=-3/10$ ($r_S=1.363$) and on Schwarzschild; with $r_S=1$ it misses by $5\times10^{-2}$ and $3\times10^{-1}$. Where $R'(\rho_h)\ne0$ the coefficient $c_2$ is pure gauge and $c_1$ is not
(\texttt{notes/black\_universe/static\_c1\_members.md}).
```

## 7. Discrepancies and changes made during the work

- #259's `static_c1.md` Sec. 8 says c₁ "was not run" at the other members, and Sec. 9 asks for this task; the run gives no disagreement with its text. Its Paper IV paragraph's "$\delta A_S=8\pi r_Sc_1$ with $r_S=1$" is the c = −1 case.
- `horizon_lemma.md`'s table of the members has no row for c = −0.3 (R′(ρ_h) = −0.679735, κ = 0.601300): proposed above.
- `full_slice.Bg(..., schw=True).kappa` is the κ of the c = −1 chart (0.6366), not 1/2; this script uses κ = 1/2. The Schwarzschild chart of `horizon_sector.jetsK` returns NaN for Π ≳ 150 (the Newton start overflows):
  the slice reads stop at s = 10.
- Bars and design changed after the first development run (all named in the raw README before the runs of record): the T7 pointwise bar 1e-3 → 5e-2 with both of the last two ratios 4.0 ± 0.3 (c = −0.3 has e^{−F} = 10.9 at S);
  the step ladder 0.01 … 0.00125 → 0.005 … 0.000625; the a(0) probe at s = 1e-6; T8's bars relative to max(|c₁|, 1e-5). The bars of the issue (1e-6, 1e-12) were not touched.
- The first chain of eight development mutations did not start (a relative path to `run_mut.sh`); a second launch produced the logs. The runs of record are one launch each.

## 8. What this does not show

- Any member but the four run. The identity is *proved* for the exact family at every Ellis member (S3), checked at c = −1, −1.5, −0.3, not at the other 17 of #212's black universes.
- That c₁ is the whole story when c₂ is gauge. At R′(ρ_h) ≠ 0 the growth of the area of the perturbed event horizon is not c₂ e^{κv} and was not computed; only the area of S is.
- Data with tails. δφ ∼ 1/ρ and a non-compact δF at infinity (the exact mode itself) are read through cuts, as in #259; the proof of δM_ADM = −κΘ needs δφ = δF = 0 beyond a finite ρ.
- Evolved data other than #248's classes and the exact mode, and the other static family members (is r = √L R the only regular static c₁ ≠ 0 perturbation of the fixed theory).
- The literature. κ δ r_S-type first laws are not checked against Bardeen–Carter–Hawking or Iyer–Wald *(unverified)*.
- Precision beyond the solver's: at c = −0.3 the evolved horizon-fit c₁ is 9.0e-7 off the exact value, 1.1 times inside 1e-6.

## 9. Recommended tasks

1. **The identity at all of #212's black universes.** *What:* `static_c1_members.py` T1 and T2 (mpmath) at the 20 members of #212, and T4 (slice code) at five of them, including one with ρ_h > 1. *Why:* three members run; the
   closed form says −c/π for all. *Acceptance:* κ r_S c₁ = −c/π and dM = κ_ξ dA/(8π) to 1e-12 at all 20; Θ_gen + r_S c₁ ≤ 1e-9 at the five; `rs_one` fails at every c ≠ −1. *Depends on:* #301, #212.
2. **The invariant area growth of the event horizon at R′(ρ_h) ≠ 0.** *What:* locate the perturbed horizon from evolved data (null-ray integration or the apparent-horizon condition of #278) at c = −1.5 and read δA(v) of its
   cross-sections against c₁, c₂ and the displacement. *Why:* c₂ is pure gauge there, so the horizon-law growth is a gauge statement; the area is the invariant. *Acceptance:* δA(v)/(8π r_S) from the horizon cross-sections is gauge
   independent to 1e-6 under the pure gauge of T6 (η(0) = 0.03) and equals c₁ at v → −∞; the displacement-free reading (c₂ alone) changes by −η(0)R′ρ′(0) under it. *Depends on:* #301, #261.
3. **The first-law reading against the literature.** *What:* read the Bardeen–Carter–Hawking and Iyer–Wald statements of δM = κδA/(8π) (and the Smarr-type relation κ r_S² + ρ_h = K found here), and record in
   `notes/literature.md` what they do and do not contain for scalar-hairy static solutions. *Why:* the reading is *unverified* in #259 and here. *Acceptance:* each equation quoted from the source read, with the match or
   mismatch in normalization of κ and M stated. *Depends on:* #260.

## 10. Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread:

```
static_c1_members.py derive [--mutate=s_q|s_pot|s_k]                                        # 3 checks, 1 s
static_c1_members.py check [--only=T1,..,T8] [--mutate=rs_one|ctl_true|c1_noshift|th_noF|th_nophi|th_asign|dm_noT|kappa_off|jump_off|xform_p|xform_q|rho_ellis|mode_khat|mode_dF0|law_norm|law_mass|gauge_noRp]   # 8 checks, 2 min
static_c1_members.py table                                                                   # the numbers of this note, 3 min
static_c1_members.py mutations DIR PREFIX                                                    # fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones
```

`static_c1_members_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I301-mac/runs`, code commit cde8be5), in that order. Memory 0.68 GB at most (the finest step of the solver, T7).

## Coordinator's review

Merged at edec471. Review folder `negative-mass-raw/I301-mac-review` at 95e1df88, criteria committed (cd63bcc6) before the first run of
`ind301.py`. The changed files are those the brief names.

- **Rerun on the Mac** from 167e537: `derive`, `check` and `table` exit 0, and their logs equal the runner's apart from timings.
- **ind301**, 4 of 4, from closed forms alone: M(L) = (K/3) L^{3/2} for every c (sympy, c symbolic) and L^{3/2}/2 for Schwarzschild;
  κ r_S c₁ = dM/dL to 2e-31 at c = −1.5, −0.3 and Schwarzschild, with the horizon's shift dρ_h/dL = 1/A′(ρ_h) in c₁; r_S = 1 misses by
  5.3e-2 and 0.267; the first law along each family to 2e-31. κ R_h² + R_h R′_h equals the static constant K at the horizon exactly,
  which is the hand proof of κ r_S c₁ = K/2. The mutations `no-shift`, `rs-one` and `no-norm` fire.
- The solver bar the runner moved at c = −0.3 (1e-3 to 5e-2) is named in its raw README before the run of record.
- Applied: `static_c1.md`, `horizon_lemma.md` (a row and a remark on c₁'s invariance), and Paper IV (`paper_numbers.py` 144 strings,
  mutation `rsone`).

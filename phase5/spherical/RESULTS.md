# Phase 5 — spherical (1+1) nonlinear evolution of a lone signed-source boson star

Model: Shirokov arXiv:2608.24577, one sector. G_ab = 8πε T_ab (ε = +1 canonical, −1 phantom), canonical
Klein–Gordon equation, V = ½|Φ|² − ¼λ|Φ|⁴ + ⅙μ|Φ|⁶, λ = 10240, μ = λ²/4.8. Static stars come from
`phase2/radial/background.py` (BVP, α(∞) = 1, true ADM masses), imported read-only.

**Short answer.** Evolved accurately in spherical symmetry, both campaign stars are stationary over t = 1500.
- Their frequency, central amplitude, Misner–Sharp mass and Noether charge stay constant. The residual drifts
  are 10⁻⁹–10⁻⁷ at h = 0.05–0.1, converge at 4th–5th order, and come from Kreiss–Oliger dissipation.
- The paper's 5%-radius seam produces only a ~3×10⁻⁵ breathing transient plus a constant ω offset of
  +3–6×10⁻⁶ (≤ 1.3×10⁻⁵ in a no-jump variant), with no drift. The result is the same at h = 0.1 and 0.05.
- Shirokov's drifts (Δω ≈ −0.02 to −0.03, content +12–21% by t = 200) are therefore not physical properties of a
  lone star, and the seam does not cause them. In this 1D setting they would have to be numerical, which I cannot
  test here: the 3D code and its gauge differ (see caveats).
- Inside these evolutions, wherever a star sits near the family, the ω→M map ("frequency drift implies mass
  change") holds to ≲5% of the drift. For a quasi-static drift it would therefore be a valid mass proxy. Whether
  the 3D drift is quasi-static is not settled by the pack's diagnostics (§6.3).

## 1. Method (`evolve.py`, `run.py`)
- **Gauge and variables.** Polar-areal gauge ds² = −α²dt² + a²dr² + r²dΩ², a⁻² = 1 − 2m/r. The evolved variables
  are the complex Φ and Π = (a/α)∂_tΦ:
  ∂_tΦ = (α/a)Π, ∂_tΠ = r⁻²∂_r(r²(α/a)∂_rΦ) − αa·2V′Φ.
- **Metric on every RK substep, from the constraints.** Only these two equations contain ε.
  - Hamiltonian constraint: m′ = 4πεr²[½K(1 − 2m/r) + V], with K = |Π|² + |Φ_r|². This m is the Misner–Sharp
    mass, which is exact and quasi-local in spherical symmetry.
  - Polar slicing: (ln α)′ = m/(r(r − 2m)) + 4πεr(½K − a²V), normalised so that α·a = 1 at R. So t is
    asymptotic proper time and ω is the physical frequency.
  - These are the standard Seidel–Suen / Hawley–Choptuik equations. The approach is fully constrained, with no
    metric evolution. That is cheap, accurate, and gives m(t, r) directly.
- **Discretisation.** Cell-centred grid r_i = (i + ½)h with even-parity ghosts at the origin; 4th-order centred
  differences; RK4 at CFL 0.25.
  - 6th-order Kreiss–Oliger dissipation (σ = 0.2, the same operator as GRTeclyn).
  - Quadratic sponge on [0.7R, R].
  - Radial ODEs integrated at 4th order: RK4/Simpson with 4-point midpoint interpolation.
- **Diagnostics** (`runs/*.csv`, every Δt = 1):
  - |Φ(0)| and its phase. ω(t) comes from linear fits to the unwrapped central phase.
  - m(r) at r = 8, 12, 20, 30, 50 and at R; the Noether charge N(r) = −∫4πr² Im(Φ̄Π) dr at the same radii.
  - Shirokov's "content" ∫√(|Φ|² + |Π|²)dV to r = 20.
  - Max momentum-constraint residual |∂_t m − 4πεr²(α/a³)Re(Π̄Φ_r)|, every 50 units.
- **Domain.** R = 300 (400 for perturbed runs). Matter radiated at t = 0 with group velocity < 1 reaches the sponge
  at t > 210. Doubling to R = 600 changes nothing (`c750_h0.1_R600` matches `c750_h0.1` to all printed digits);
  so does halving to R = 150 at h = 0.5.
- **Cost.** h = 0.1, R = 300, t = 1500 takes ~3 min on one core (numba).

## 2. Validation

### 2a. Unstable-branch growth rates (`growth.py`, `growth_fit.py`, `growth_output.txt`)
Setup: ω = 0.94 stars seeded with Φ, Π × (1 ± δe^{−r²/4}). I fit the odd part D(t) = (|Φ_c|₊ − |Φ_c|₋)/(2|Φ_c|₀)
to Ae^{σt} + Be^{−σt} + C. Using the odd part cancels the grid-seeded perturbation and O(δ²) terms.

| star | radial σ = √(−σ²) (phase2) | h = 0.1, fit [0,80] | h = 0.05, fit [0,80] | δ = 1e-7 vs 1e-5 |
|---|---|---|---|---|
| canonical ω = 0.94 | 0.07370 | 0.07332 (−0.5%) | **0.07376 (+0.08%)** | identical to 5 digits |
| phantom ω = 0.94 | 0.06668 | 0.06635 (−0.5%) | **0.06679 (+0.16%)** | identical to 5 digits |

On window [0,60] the fits come out at +0.2–0.5%. The residual depends on the window at the 0.5% level because of
the slowly decaying continuum content.

**The nonlinear code reproduces the linear radial growth rate for both signs of ε.** That independently confirms
the phase2 radial eigen-solver, including how ε enters.

Nonlinear end states (h = 0.1):
- δ > 0 migrates to a stable star with ω ≈ 0.92 (canonical) / 0.914 (phantom), |Φ_c| up ~10–13%, and m(R)
  conserved to 1e-5.
- δ < 0 disperses: |Φ_c| falls 99%, ω → 1, m(20) → 4% of M.

### 2b. Convergence of the stable campaign stars (`drift_table.py` → `drift_table_output.txt`, `drift_table.csv`)
Exact BVP seed, R = 300, t = 1500. The drift is measured between windows [0,100] and [1400,1500].

| run | h | ω[0,100] − ω_star | Δω | Δ\|Φ_c\|/\|Φ_c\| | Δm(50)/m | ΔN(50)/N | m(50)/M_ADM − 1 |
|---|---|---|---|---|---|---|---|
| canonical 0.750 | 0.5 | −1.42e-04 | +2.20e-04 | −9.1e-05 | −1.16e-03 | −1.42e-03 | −4.7e-04 |
| | 0.333 | −3.06e-05 | +3.29e-05 | −1.5e-05 | −1.74e-04 | −2.13e-04 | −9.5e-05 |
| | 0.25 | −9.88e-06 | +9.11e-06 | −4.5e-06 | −4.83e-05 | −5.92e-05 | −3.0e-05 |
| | 0.2 | −4.10e-06 | +3.27e-06 | −1.6e-06 | −1.74e-05 | −2.13e-05 | −1.2e-05 |
| | 0.1 | −2.62e-07 | +1.17e-07 | −6.6e-08 | −6.29e-07 | −7.71e-07 | −7.7e-07 |
| | 0.05 | −1.66e-08 | +3.75e-09 | −2.6e-09 | −2.05e-08 | −2.51e-08 | −4.8e-08 |
| phantom 0.7603 | 0.5 | −1.59e-04 | −3.31e-05 | +1.3e-05 | +1.72e-04 | +2.11e-04 | −4.0e-04 |
| | 0.333 | −3.36e-05 | −4.58e-06 | −3.4e-07 | +2.36e-05 | +2.88e-05 | −8.3e-05 |
| | 0.25 | −1.08e-05 | −1.57e-07 | −6.3e-07 | +4.8e-07 | +6.0e-07 | −2.7e-05 |
| | 0.2 | −4.44e-06 | +1.74e-07 | −3.3e-07 | −1.10e-06 | −1.33e-06 | −1.1e-05 |
| | 0.1 | −2.81e-07 | +1.74e-08 | −2.2e-08 | −1.06e-07 | −1.29e-07 | −7.1e-07 |
| | 0.05 | −1.76e-08 | +5.9e-10 | −1.2e-09 | −4.0e-09 | −4.9e-09 | −4.4e-08 |

- **Observed orders.** The ω offset and the m(50)/M_ADM error converge at order 3.97–4.01 at every step, as
  designed. The canonical drift Δω converges at order 4.5–5.0. The phantom drift passes through zero near
  h ≈ 0.25, which is where the odd orders in the table come from. Its finest pair is order 4.9.
- **Momentum constraint.** The max residual is 1.1e-7 / 4.0e-9 / 1.3e-10 at h = 0.2 / 0.1 / 0.05 (order 4.8–5),
  and it is flat in time on every run.
- **What drives the drift: KO dissipation, not the boundary.**
  - At h = 0.5, setting σ_KO = 0 reduces Δω from +2.2e-4 to +1.6e-5 (canonical) and from −3.3e-5 to +1.7e-5
    (phantom).
  - At h = 0.5, σ_KO = 1 raises it to +9.1e-4 / −1.7e-4.
  - At h = 0.1, σ_KO = 0 gives Δω = +3e-9 against 1.2e-7 with σ_KO = 0.2 (canonical), and +3e-9 against 1.7e-8
    (phantom). Δm(50) without KO is −3e-8 for both.
  - The residual drift is a slow, dissipation-driven loss of N, which moves the star along its family. For ε = +1
    and h = 0.5 it is a loss; ε = −1 shows a small gain. The sign depends on ε, but both are ≲10⁻³ in N.

### 2c. Mass and charge diagnostics are consistent
- **Misner–Sharp mass vs ADM mass.** At h = 0.05 and t = 0, m(r)/M_ADM − 1 is
  −4.5e-4 at r = 8, −2.3e-6 at r = 12, and −4.7e-8 at r = 20, 30, 50 and R (the matter tail).
  Over t = 0–1500 its maximum deviation is 2.2e-8 (canonical) / 4.3e-9 (phantom). See `misner_sharp_output.txt`.
- **Noether charge.** N(50) reproduces the BVP N to 4.5e-10.
- **Phantom sign.** For ε = −1 the Misner–Sharp mass is negative everywhere outside the core and conserved like the
  canonical one.

## 3. Lone campaign stars to t = 1500 (question a)
Both stars are stationary to truncation error. At h = 0.1 (5× finer than Shirokov's finest dx = 0.25), over
t = 0 → 1500:

| | canonical ω = 0.750 | phantom ω = 0.7603 |
|---|---|---|
| \|Φ_c\| drift | −6.6e-8 | −2.2e-8 |
| ω drift | +1.2e-7 | +1.7e-8 |
| Misner–Sharp m(50) drift | −6.3e-7 | −1.1e-7 |
| N drift | −7.7e-7 | −1.3e-7 |

Even at h = 0.5 (Shirokov's coarsest dx) all drifts over t = 1500 are ≤ 1.4e-3. That is 4–5 orders of magnitude
below Shirokov's numbers, all of which are at t = 200 and at dx ≤ 0.5:
- phantom peak −0.63% per 100 time units;
- content +12–21%;
- canonical ω 0.748 → 0.728.

Comparison of time series: `mass_vs_time.csv` (m at 8, 12, 20, 30, 50 and R; N; ω(t); |Φ_c|).

Perturbed stable stars (`endstate_check.py` → `endstate_check_output.txt`, `proxy_test_output.txt`):
- **5% central bump** (Φ, Π × (1 + 0.05e^{−r²/4})): both stars radiate 0.13% of M and 0.04% of N within
  t ≈ 150, then sit on a neighbouring static star. The remnant at the measured late ω matches the BVP star at that
  ω to 5e-7 in |Φ_c| and 2e-6 in m and N. After t ≈ 300 nothing drifts.
- **Uniform ±3% amplitude rescale:** same behaviour. The star relaxes within ~100 units to ω shifted by ∓0.009
  and stays there; late breathing is ≤ 4e-7.
- Phantom and canonical respond almost identically (Q-ball-dominated regime, |M|/R ~ 3e-3).

## 4. Misner–Sharp mass vs ω: test of the frequency→mass proxy (question b)
Family slopes come from the BVP neighbours (`family_slopes.py`):
- canonical at 0.750: dlnM/dω = −5.28, dlnN/dω = −6.47;
- phantom at 0.7603: dlnM/dω = −5.34, dlnN/dω = −6.51;
- first law dM = εω dN holds to 4e-4.

I compared the predicted ΔM/M = (dlnM/dω)·Δω with the measured Δm(50)/m in every run where the drift is resolved:

| run | Δω | predicted ΔM/M | measured Δm(50)/m | ratio |
|---|---|---|---|---|
| canonical h = 0.5 | +2.20e-4 | −1.16e-3 | −1.16e-3 | 1.002 |
| canonical h = 0.5, σ_KO = 1 | +9.10e-4 | −4.81e-3 | −4.79e-3 | 0.995 |
| canonical h = 0.333 | +3.29e-5 | −1.74e-4 | −1.74e-4 | 1.003 |
| canonical h = 0.2 | +3.27e-6 | −1.73e-5 | −1.74e-5 | 1.006 |
| phantom h = 0.5 | −3.31e-5 | +1.77e-4 | +1.72e-4 | 0.971 |
| phantom h = 0.5, σ_KO = 1 | −1.74e-4 | +9.30e-4 | +9.25e-4 | 0.995 |
| canonical seam h = 0.1 (t ≥ 200) | +1.0e-7 | −5.4e-7 | −5.4e-7 | 1.006 |

The ratio departs further from 1 (0.88–1.5) only where Δm ≲ 2e-6, i.e. at the level of residual breathing and
truncation.

The proxy also holds for relaxed perturbed stars, comparing the shift from the original star with the static
family:
- +3% rescale: measured m/M − 1 = +4.92% against +4.74% from (dlnM/dω)·Δω. The linear estimate's error at
  Δω = −0.009 is the curvature of M(ω).
- 5% bump: +2.89% against +2.83%.

**Conclusion on the proxy.** In spherical symmetry, a star that drifts slowly and stays close to the static family
(numerical dissipation is one example) carries the gravitating mass M(ω(t)) to within a few percent of the drift.
So the proxy used in `review/acceleration_drift` is internally sound *for a quasi-static drift*.

Limits of the test:
- It only covers drifts 10²–10⁴× smaller than the 3D one.
- It cannot tell whether the 3D stars are close to the family (§6.3).

## 5. Seam test (question c-ii) (`seam_analysis.py`)
Seed as in GRTeclyn `boson_star_ode.py`:
- the BVP core up to r_m (first φ < 5%φ_c: r_m = 5.19 / 5.22);
- beyond r_m, the analytic tail φ_m(r_m/r)e^{−κ(r−r_m)} and a Schwarzschild lapse with the enclosed mass at r_m
  (0.0143454 / −0.0142938);
- the core lapse rescaled so that α jumps by −4.3e-5 / +5.1e-5 at r_m (F9 values);
- the lapse enters only through Π = −iω(a/α)φ, and the metric is re-solved from the constraints, as in GRTresna.

| | canonical (h = 0.1) | phantom (h = 0.1) |
|---|---|---|
| constraint-solved M vs exact seed | −3.1e-5 (rel) | −1.8e-5 |
| N vs exact seed | −3.8e-5 | −2.2e-5 |
| ω (t = 500–1500) − exact-seed ω | +5.82e-6 | +3.41e-6 |
| predicted from ΔN via dN/dω | +5.81e-6 | +3.40e-6 |
| \|Φ_c\| breathing, t < 100 (exact seed / seam) | 1.7e-6 / 2.8e-5 | 1.9e-6 / 2.7e-5 |
| breathing at t = 1400–1500 | 1.7e-9 / 2.0e-8 | 3e-10 / 7e-9 |
| drift [0,100] → [1400,1500] in ω, \|Φ_c\|, m, N | ≤ 8e-7 | ≤ 2e-7 |

What the seam does:
- It seeds a star with 2–4e-5 less charge.
- That star settles immediately onto the static neighbour with that charge. The ω offset matches the family
  prediction to 0.3%.
- It adds a small decaying breathing mode.

It produces no relaxation or drift resembling the paper's, so on its own it does not explain Shirokov's drifts. The
seam-mass/ADM-mass mismatch (1.9%) does not show up in the evolved star's mass, because the constraint solve fixes
the metric: the evolved star carries the true ADM mass to 3e-5. This matches the F9 expectation.

**The seam effect is resolution-independent (it is in the data, not in truncation).** At h = 0.05 both seams
repeat the h = 0.1 numbers:

| h = 0.05 | ΔN | Δω (predicted) | early breathing | drift ω / m |
|---|---|---|---|---|
| canonical | −3.75e-5 | +5.82e-6 (+5.80e-6) | 2.8e-5 | +9e-8 / −6e-8 |
| phantom | −2.21e-5 | +3.41e-6 (+3.40e-6) | 2.6e-5 | +9e-8 / −4e-8 |

The mean |Φ_c| offset settles by t ≈ 100 and is then flat to 1e-8 through t = 1500, in windows
(−3.3, −2.8, −2.8, −2.8, −2.9)e-6 at h = 0.1.

**No-jump variant.** Same tail and seam mass, but the core lapse is rescaled to be continuous at r_m. This rescales
α, and therefore Π, everywhere in the core, so it is a different seed, not "the tail alone".
- canonical: ω − ω_exact = +1.24e-5 (with jump +5.8e-6);
- phantom: −4.3e-6 (with jump +3.4e-6);
- breathing 2.8e-5 / 3.4e-5.

So the size of the offset depends on how the core lapse is normalised. In every variant it is ≤ 1.3e-5 in ω, more
than 3 orders of magnitude below the paper's drift, and constant in time.

## 6. What this says about Shirokov's drifts
1. **Not physical.** A lone stable star, evolved accurately with the full nonlinear (spherical) Einstein–KG system
   for either sign of ε, is stationary to ≤1e-6 over t = 1500. Seeding errors much larger than the paper's seam
   (5% bumps, 3% rescales) relax within ~150 units onto a nearby static star and then stop. A lone star has no
   physical mechanism for a monotone ω drift of 0.02–0.03 by t = 200.
2. **Not the seam.** See §5.
3. **Is the 3D drift a quasi-static slide along the family? Inconclusive** (`pack_compare.py`, `family_diag.py`).
   If the pack's stars slid along the static canonical family from ω_Q1 to ω_Q4, the static stars predict:

| grid | ω (Q1→Q4 midpoints, ~147 units apart) | static ΔM/M | static Δcontent | reported Δcontent (t 0→200) | static Δpeak | reported Δpeak |
|---|---|---|---|---|---|---|
| N128 | 0.7451→0.7148 | +18.6% | +13.7% | +21.0% | −0.60% | −0.36% |
| N192 | 0.7466→0.7218 | +14.8% | +10.7% | +16.1% | −0.44% | −0.34% |
| N256 | 0.7479→0.7283 | +11.3% | +8.1% | +12.4% | −0.30% | −0.25% |
| lone N128 | 0.7433→0.7135 | +18.3% | +13.6% | +20.0% | −0.62% | −0.48% |

   - **Content.** The reported content change is 1.5× the static prediction on every grid. However, the content
     spans t = 0→200 and the ω difference spans about 147 units. For a linear drift that window mismatch alone gives
     a factor of ≈ 1.36, and the pack's content definition (Cartesian coordinate volume, its Π normalisation) is
     only reconstructed here. So the content numbers are **roughly consistent** with a slide.
   - **Peak.** It changes by 60–85% of the static prediction. This is also not decisive, since the window of the
     peak numbers is uncertain.
   - **Metric.** Against the slide: the N128 lone star's 1−χ at R = 8/16 grows only 1–2% over t = 40–200
     (review/acceleration_drift H1), while a slide needs ≈ +15% over that span. χ is gauge-dependent, so this is
     indicative only.
   - **Net.** The pack diagnostics do not establish whether the 3D stars stay on the static family.
4. **So the drifts are most plausibly numerical artefacts of the 3D evolution** (or of its initial-data solve).
   They are absent from an accurate evolution of the same physics in spherical symmetry. Plausible candidates, all untested here:
   - KO dissipation and sponge at dx = 0.25–0.5. In 1D, σ_KO acts as a mass-loss channel, but even σ_KO = 1 at
     h = 0.5 gives only 5e-3 by t = 1500.
   - Under-resolution of the steep core (φ falls by e from r ≈ 2.5 to 3.5), which is worse in 3D Cartesian at
     dx = 0.5 than in 1D at h = 0.5.
   - The CCZ4/moving-puncture gauge and constraint damping.
   - The non-convergent order 0.5–0.9 the review found is itself the sign of an error that is not in the
     asymptotic regime.
5. **Implication for the frequency→mass proxy.**
   - The proxy is correct physics for quasi-static changes. §4 verifies it to ≲5% of the drift in the
     nonlinear 1D evolution, including relaxed stars 3–5% away in mass.
   - Whether it applies to the 3D runs depends on whether the 3D star stays near the family. Point 3 leaves that
     open. The N128 1−χ data hint that it does not.
   - So the review's "unconfirmed" stands. A per-star metric mass vs t in 3D would settle the fine-grid
     acceleration rise: an ADM-like surface integral or the areal mass on spheres around a lone star, on the
     N192/N256 grids.

## 7. Caveats
- Spherical symmetry only. It says nothing about l ≥ 1 dynamics, the pair, or 3D-specific numerical modes.
  Radial stability plus 1D nonlinear stationarity establish that the *physics* has no drift, not where the 3D
  code's drift comes from.
- Polar-areal gauge versus Shirokov's 1+log/Γ-driver: a gauge cannot create a monotone drift in a
  gauge-invariant ω measured at infinity, but the 3D numbers mix coordinate time and coordinate volume.
- The frequency is measured from the phase of Φ at r = 0 in asymptotic time (α·a = 1 at R). In polar slicing α(0)
  is constant for a stationary star, so this is the physical ω.
- No AMR or Cartesian-grid effects are reproduced. My h = 0.5 is not equivalent to dx = 0.5 in 3D: Cartesian
  stencils resolve the radial profile anisotropically.
- The comparison in §6.3 uses the pack's diagnostics as transcribed in review/acceleration_drift and my
  reconstruction of their definitions (modulus √(φ₁² + Π₁² + φ₂² + Π₂²), Π = ω|Φ|/α for the static star, areal
  volume). The time windows of ω and content differ. Treat §6.3 as order-of-magnitude.
- The phantom drift at coarse h changes sign with resolution, so its coarse-grid sign has no meaning.

## Files
| file | role |
|---|---|
| `make_stars.py` → `stars/*.npz`, `make_stars_output.txt` | BVP stars (campaign, ω = 0.94, family neighbours) |
| `evolve.py` | 1+1 solver (numba) |
| `run.py`, `batch.sh`, `jobs_*.txt` | runs with checkpoint/resume → `runs/<name>.csv` (+ `.final.npz` profiles) |
| `analyze.py`, `drift_table.py` → `drift_table_output.txt`, `drift_table.csv` | drift and convergence tables |
| `growth.py`, `growth_fit.py` → `growth_output.txt` | unstable growth rates |
| `family_slopes.py` → `family_slopes.json/_output.txt` | dM/dω, dN/dω |
| `proxy_test.py` → `proxy_test_output.txt`, `mass_vs_time.csv` | ω→M proxy test; Misner–Sharp m, N, ω(t) time series |
| `endstate_check.py` → `endstate_check_output.txt` | relaxed perturbed stars vs static family |
| `seam_analysis.py` → `seam_analysis_output.txt` | seam test |
| `family_diag.py`, `pack_compare.py` → `family_diag_output.txt`, `pack_compare_output.txt` | pack drifts vs quasi-static slide |
| `misner_sharp_output.txt` | m(r) at several radii vs M_ADM |

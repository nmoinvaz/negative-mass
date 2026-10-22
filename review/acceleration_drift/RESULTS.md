# Acceleration drift on the fine d=10 grids — RESULTS

Data: GRTeclyn pack (commit 09fbcd6), `results/bondi-dipole-runaway/campaign`. I only read the data and executed nothing from
the repo. Scripts in this directory, run as `python3 -I <script> <campaign_dir> .`; their stdout is in `*_output.txt`.
Masses: M+ = 0.0146265, |M−| = 0.0145803, M̄_true = 0.0146034. "a" is a quadratic fit to the core-tracker midpoint
(`sector_dynamics.dat` cols 2, 6) unless stated otherwise. I reproduced F2's quarters exactly (N128 1.468/1.445/1.460/1.469,
N192 1.478/1.524/1.594/1.662, N256 1.480/1.522/1.580/1.641 ×1e-4).

## What the pack offers as a "mass" diagnostic
The pack has no per-sector ADM, Komar or Noether-charge diagnostic. The candidates are:
- `sector_barycenters.dat` total_canon/total_phantom = ∫ sqrt(φ₁²+Π₁²+φ₂²+Π₂²) dV (coordinate volume; this is the
  `_matter_sectors` modulus). It is linear in amplitude and is not a charge or an energy. `confinement.dat` gives the fraction
  within r = 6 of each sector's barycentre.
- `sector_dynamics.dat` peak_canon/peak_phantom: the maximum of the same modulus.
- `collapse_diagnostics.dat` max_phi + min_phi (and max_Pi + min_Pi): this is the global extreme of the canonical real part.
  Its zero crossings give the star's **coordinate oscillation frequency ω(t)**. It is **canonical sector only** (consistent
  with cpp F6). I turned ω(t) into a mass proxy through the equilibrium family M(ω): dM/dω = −0.0661 near ω = 0.75
  (phase2 `boson_star_output.txt`). Justification: for a star that stays near equilibrium, ω labels its position on the
  family. Caveat: this is coordinate time, but min lapse changes by only 0.06% over t = 5–200, so lapse redshift does not
  explain the shift.
- `shell_profiles.dat`: the 1−χ shell mean at R = 8 and 16. This is a metric (gravitating) mass proxy, but it is usable only
  for the box-centred lone canonical star. For the pair the dipole cancels it.
- Matter momentum p/v (from `px_*`): an effective inertial mass.

## H1 — mass growth (verdict: **inconclusive; the field-content version is ruled out**)
| grid | content growth t=0→200 (c / p) | peak | ω, Q1→Q4 | M_ω/M_ω(Q1) at Q4 | a_Q4/a_Q1 | (a/M_ω) Q4/Q1 |
|---|---|---|---|---|---|---|
| N128 | +21.0% / +21.4% | −0.36% | 0.7451→0.7148 | 1.138 | 1.000 | **0.879** |
| N192 | +16.1% / +15.8% | −0.34% | 0.7466→0.7218 | 1.112 | 1.125 | 1.012 |
| N256 | +12.4% / +12.0% | −0.25% | 0.7479→0.7283 | 1.088 | 1.109 | 1.019 |
| lone canonical N128 | +20.0% | −0.48% | 0.7433→0.7135 | — | — | — |

- Field content grows **most on N128, where a stays flat**, and grows less on the finer grids. Across grids the content
  growth runs opposite to the acceleration rise, so content growth cannot be the cause. The lone N128 stars grow by the same
  20–21%, which shows it is a single-star evolution effect and not something the pair causes. Content and ω-drift converge at
  order only ≈0.5–0.9.
- With the ω-based proxy, a/M stays constant to within 1–2% on both fine grids. Ranking across the fine pair also agrees:
  N192 has more ω-drift and a larger rise than N256. **On N128 the same proxy fails by 12%.** In the late N128 cells
  (`late_ad2.csv`), a·d²/M_ω sits at 0.88–0.92 for t ≈ 150–600 in chase03c, while d goes 10→17 and a falls by 2.3×. That
  is suggestive but includes the early 12% failure.
- Against the metric proxy, the lone canonical N128 star's 1−χ grows only +1.2% at R16 and +2.0% at R8 over t = 40–200,
  versus the ≈+14% that ω implies. So on N128 at least, the ω drift and content growth do **not** show up as a comparably
  larger gravitating mass. χ is gauge-dependent, so this is indicative rather than conclusive. The pack has no fine-grid lone
  star, so the metric check cannot be repeated on N192/N256.
- Conclusion: a "growing mass" explains the fine-grid rise numerically only with the ω proxy. The same proxy is contradicted
  by N128, and N128's metric does not show the growth. It is not established.

## H2 — gauge / kinematics (verdict: **ruled out**, total ≲0.3% against an 11% rise)
- v_mid(200) = 0.028–0.030, so γ³−1 = 0.12–0.14% at the end and smaller averaged over Q4.
- The mean shift at the cores is ≈ −1×10⁻⁵. Its drift rate is −1×10⁻⁷ to −2.4×10⁻⁷, i.e. 0.07–0.16% of a. The coordinate
  grid is not drifting.
- Min lapse: 0.9918 → 0.9913 on every grid (0.06%).
- The 1PN midpoint has no velocity correction to proper acceleration (phase3).
- The gauge twin (η = 2, N128) reproduces the quarters to 0.1%. It exists only at N128.

## H3 — separation change (verdict: **ruled out**)
On the fine grids the core separation stays within 9.997–10.014 and the proper separation within 10.000–10.024 over all
quarters. Using the instantaneous d changes a·d² by ≤0.5%. a·d²/M̄_true by quarter:

| grid | fixed d=10 | core d | proper d |
|---|---|---|---|
| N128 | 1.005/0.990/1.000/1.006 | 1.005/0.989/1.002/1.019 | 1.006/0.990/1.004/1.017 |
| N192 | 1.012/1.044/1.091/1.138 | 1.012/1.043/1.091/1.142 | 1.012/1.044/1.093/1.144 |
| N256 | 1.014/1.042/1.082/1.124 | 1.014/1.042/1.082/1.125 | 1.014/1.042/1.083/1.127 |

The estimator disagreement in cpp F4 is ≈0.1 at most, which is far too small to matter.

## H4 — numerical (verdict: **the rise is not a converging error at the available resolutions; sponge/box and constraint growth are ruled out at N128 and not testable on fine grids**)
- **Convergence.** Q4/Q1 is 1.0004 (N128), 1.1250 (N192), 1.1086 (N256). The fine pair differs by only 1.6 points of an
  11-point rise. If the rise were truncation error it would need order p ≈ 0.5 to vanish. Richardson at p = 2 keeps the rise
  (Q4/Q1 → 1.088). N128 is the outlier, not the fine grids.
- **Symmetric vs antisymmetric** (`sym_anti.csv`). N128 develops a large antisymmetric acceleration (canonical +10.5%,
  phantom −10.4% Q1→Q4; anti/sym reaching 0.10). This is the pack's "opening" artefact, which converges at order ≈3. On the
  fine grids it is 0.035/0.017 of sym. Its symmetric part stays flat. Notably, N128's canonical star alone rises like the
  fine-grid midpoint. One possibility is that N128's flatness comes from two compensating errors. This is speculation; the
  data cannot separate them.
- **Constraints.** L2 Ham changes ×1.01–1.08 between t = 5 and 200 on all L64 grids, and Mom falls. There is no growth.
- **Sponge / box.**
  - At t = 200 the leading (canonical) core is at r = 8.0 from the box centre, i.e. 16 from the sponge's inner radius (24).
    The φ = 5%·φ_c radius (≈5.2) puts the star's visible edge ≈11 from the sponge. Net boundary flux is ~1e-11 (round-off)
    after t ≈ 120 on all grids.
  - The doubled box (L128, N256, same dx = 0.5 as N128) matches N128's quarters to 0.2–0.4%.
  - nomill_left (pair offset 10 to the left) matches the centred longrun's drift to 0.1% at t = 200 and 0.5% at t = 400.
  - So the box/sponge does not affect the N128 acceleration. All of these are dx = 0.5 runs, so they say nothing about the
    fine grids, where the sponge damping is ~30× weaker (cpp F3). No fine-grid box test exists.
- Elliptic solve depth (deepsolve) and one AMR level (amrcheck) are identical to N128 to <0.1%.

## Late-time cells (N128 only)
- **longrun t=400.** The midpoint a stays flat (1.45–1.48e-4) only to t ≈ 300, then falls: 1.379 (305–355) and 1.264
  (350–400). Its 1.418e-4 "last ⅔" fit averages this out. Meanwhile d opens 10→11.8, and the antisymmetric part grows until
  the phantom's own acceleration is 0.75e-4 against the canonical's 1.78e-4.
- **nomill_left (t=600)** and **chase03c treadmill (t=784)** reproduce the same curve: a falls to 0.64e-4 at d = 15 and
  0.22e-4 at d = 26. a·d²/M̄_true rises to 1.16–1.31, so with fixed mass a does not follow 1/d² during the opening either.
- "Steady to 2%" holds at N128 only for t ≲ 300, and only for the midpoint. It also holds only because d stayed ≈10.

## Best explanation (hedged)
The fine-grid rise is not caused by:
- gauge or relativistic kinematics (≤0.3%),
- separation change (≤0.5%),
- constraint growth,
- the box at dx = 0.5.

It barely changes between N192 and N256, so at these resolutions it behaves like a resolved feature of the simulation, not
like convergent truncation error. The only diagnostic that tracks it quantitatively is the stars' slow downward drift in
frequency (ω 0.748→0.728 on N256). In equilibrium terms that means each star climbs the M(ω) branch by ≈9–11%, and a/M_ω is
then flat to 2% on both fine grids. However, this drift is a non-equilibrium single-star effect that also appears in the lone
controls. It is largest on N128, where it does **not** raise the acceleration, and there the metric does not register a
matching mass increase. My best reading:
- the stars are not stationary on the grid (the slow drift in ω and content is a numerical, resolution-dependent evolution of
  the stars);
- on the fine grids the gravitating mass most likely does grow along with it;
- on N128 a large, order-3 antisymmetric artefact (the opening) plausibly masks the same effect in the midpoint.

This is consistent with the data but not proven. A fine-grid lone star with 1−χ or ADM-surface mass vs t would settle it.

## Implications for the paper
- **"a = M̄/d², steady to 2%"** is an N128-only statement. On both fine grids a rises 11% in t ∈ [5,200] at fixed d. At
  N128 it holds only to t ≈ 300 (longrun, nomill, chase all decline afterwards).
- **The 6–8% fine-grid excess** (converged d=10 a·d²/M̄_true ≈ 1.06) is the **time average of a rising quantity**, not a
  steady value:
  - In the first quarter (t = 5–54) the fine grids give a·d²/M̄_true = 1.012 (N192) and 1.014 (N256); Richardson p=2 gives
    ≈1.017.
  - By Q4 it is 1.12–1.14.
  - If the rise is mass growth, the early value (≈1.4% above the initial true-mass law, comparable to the 1PN ±0.14% body
    terms plus the ~1% scatter in fit windows) is the relevant test of a = M/d², and the "excess" reflects the stars not
    being stationary.
  - If it is not mass growth, then a(t) at fixed d is not constant and the force law is not demonstrated on the fine grids.
- Either way, the headline 1.549×10⁻⁴ / "8% high" depends on the fit window, and the paper should report a(t) and a
  per-star mass diagnostic vs t on the fine grids.

## Files
- Scripts: `drift_analysis.py`, `freq_late.py`, `h1_omega_mass.py`, `late_ad2.py`, `sym_anti_chi.py`, `per_star.py`,
  `probe.py`.
- Outputs: the matching `*_output.txt` files.
- CSV for plots:
  - `quarters.csv` — per quarter: a, per-star a, d, a·d²/M, content, peak, p/v
  - `a_mid_sliding.csv` — sliding 40-unit a(t), all d=10 cells
  - `per_star_sliding.csv`
  - `omega_windows.csv`
  - `h1_omega_mass.csv`
  - `late_windows.csv`
  - `late_ad2.csv`
  - `sym_anti.csv`

## Caveats
- The ω proxy rests on canonical-sector data only, and on an equilibrium M(ω) slope from our own solver applied to a star
  that is not in equilibrium.
- The quadratic fits on 49-unit quarters carry ~±1% window noise. The paper's barycentre windows are noisier still (F6).
- I did not inspect the frame PNGs.

# negative-mass

Working notes and code on negative (active) mass in general relativity: Bondi's 1957 runaway, the
stability of negative-mass Schwarzschild, regular negative-mass sources, and what remains open.

Repository: https://github.com/nmoinvaz/negative-mass

## Paper I — point-mass thresholds
N. Moinvaziri, *The instability of negative-mass Schwarzschild is a point-mass effect: exact thresholds for finite
stars* — [`paper/1-thresholds/main.pdf`](paper/1-thresholds/main.pdf)
- **Exactly one unstable mode.** The Dotti–Gleiser operator has exactly one negative eigenvalue for every l ≥ 2, at
  their growth rate k/(2|M|) (Theorem 1).
- **Exact count for a finite star.** With a Robin surface condition κ, the number of unstable modes is
  [R < r_s] + [κ < κ_c(R)], with κ_c in closed form (Theorem 2).
- **Thresholds.** A star with |M|/R < 0.439 and κ ≥ 0 has no unstable even-parity mode for any l; any finite radius
  caps the l⁴ growth of the point mass.
- **Verified by:** `phase2/dg_proof/` (exact SymPy checks of every identity, symbolic in l and for l = 2–10),
  `phase2/dg_spectrum/` (independent numerical spectra, 1000-case scan).

## Paper II — negative-mass boson star
N. Moinvaziri, *Stability of a negative-mass boson star* —
[`paper/2-boson-star/main.pdf`](paper/2-boson-star/main.pdf)
- **Radially stable.** Shirokov's phantom star (ω = 0.7603) is stable; the family loses stability where the Noether
  charge N peaks. The first law dM = −ω dN makes the criterion flip sign when written in M.
- **Nonlinear evolution.** 1+1 simulations keep stable stars stationary to t = 1500 and reproduce linear growth rates
  of unstable ones to < 0.2 %.
- **Non-radial.** No unstable l = 2, 3 polar or axial mode; the star supplies κ_eff > κ_c, so the point-mass
  instability is not triggered. Numerical scans, not proofs.
- **Verified by:** `phase2/radial/`, `phase2/nonradial/`, `phase2/mass_check/`, `phase5/spherical/`.

## Paper III — ghost-scalar wormholes at large D
N. Moinvaziri, *The instability of ghost-scalar wormholes at large dimension* —
[`paper/3-wormhole-largeD/main.pdf`](paper/3-wormhole-largeD/main.pdf)
- **Growth rate.** The Ellis–Bronnikov wormhole's single unstable radial mode has
  W²a² = 2(D − 3) − 4 ln 2 + (π²/3 + 4 ln²2)/(D − 3) + O(D⁻²), by matched asymptotics (the third term's closed form
  identified numerically).
- **Mode.** At large D the throat moves rigidly and the mode is (r/a) e^{1 − r²/a²} in the nearly flat region outside.
- **Threshold.** 1/W is the slope of the type I lifetime law between collapse and expansion; published runs in D = 4, 5, 6
  are consistent with it.
- **Verified by:** `notes/wormhole_largeD/` (symbolic check with control, two integrations against Torii–Shinkai's table,
  `paper_numbers.py`).

All three papers: AI-assisted (see each paper's "How the results were obtained"); not yet reviewed by specialists.

## Repository layout
- `PLAN.md` — phases, status, next steps
- `paper/` — the two papers above, with `figures/make_figures.py` regenerating every figure from repo data
- `notes/literature.md` — papers read, with what each establishes
- `phase0/` — verification of published perturbation results for negative-mass Schwarzschild
- `phase2/` — regular negative-mass (phantom boson) stars: backgrounds, Dotti–Gleiser spectrum, radial and
  non-radial stability, mass-convention check against Shirokov
- `phase3/` — first post-Newtonian dynamics of the Bondi dipole
- `phase4/` — Witten-type lower bounds on ADM/Bondi mass with a signed source
- `review/` — code review of Shirokov's GRTeclyn (arXiv:2608.24577): `python/` (initial data, analysis),
  `cpp/` (evolution)

## Findings so far

1. **Negative-mass Schwarzschild is linearly unstable** (Dotti–Gleiser 2006/2008, reproduced). The earlier
   stability claim (Gibbons–Hartnoll–Ishibashi 2004) rests on a perturbation whose metric is discontinuous.
   Exactly one unstable mode per l — now **proven for all l** (phase2/dg_proof), localized at the singularity, growth ∝ l⁴.
2. **A regular star hides that instability.** φ₀ = 1/ψ₀ is an exact zero mode of the transformed operator, giving
   a closed-form threshold: with a non-energy-injecting surface condition (κ ≥ 0) a star is stable once
   R > r_*(l) (l = 2: 2.279|M|). The l⁴ growth belongs to the point singularity only.
3. **Shirokov's phantom boson stars reproduced.** Profiles match his tables to 1.4×10⁻⁵; his "ADM masses" are
   the mass enclosed where φ = 5 % φ_c and run 1.5–2.3 % below the true ADM mass.
4. **The phantom campaign star is radially stable.** The family loses stability exactly at dN/dω = 0; the
   turning-point criterion carries over to negative mass when stated in N (not M).
5. **The phantom point particle is Bondi's all-negative mass**: one parameter is inertial, passive and active
   mass; gravity alone cannot distinguish "positive inertial, negative active". At 1PN the dipole midpoint
   acceleration is exactly m/d² (no correction); a matched pair slowly contracts (d̈ = −E_tot/d²) — the separation law and rigidity condition were
   already found by Ivanov 2012 (arXiv:1204.1269).
6. **Mass bounds:** per slice, M ≥ −(phantom's weighted energy), tight to 0.1 % for the phantom star; no bound
   uniform in time. Radiation (phase4/radiation): the hyperbolic ± pair radiates at a constant rest-frame rate
   L₀ = (24/35)m⁶/d⁶ (checked two ways), so eternal motion loses unbounded energy, but in the exact
   Bonnor–Swaminarayan solution it arrives only at the endpoint of null infinity; a 1PN back-reaction toy model
   suggests a pair from rest stops accelerating, giving a finite (huge) loss. M_B → −∞ is not established.
7. **Review of Shirokov's code:** the runaway physics is implemented correctly (signs, couplings, diagnostics
   verified). Two real bugs off the paper's path (boosted-star Π sign; Ψ4 matter term zeroed). Precision
   claims are weaker than stated: the a = M̄/d² agreement rests on the coarse grid (fine grid 7 % higher at
   d = 10), and the acceleration rises ~11 % over the run on both fine grids. That rise is not gauge, separation,
   constraint or boundary error and barely converges between the two fine grids; field-content growth runs
   the wrong way across grids. Only the stars' slowly drifting frequency tracks it (implying ~9–11 % mass growth),
   unconfirmed. Early-time fine-grid a·d²/M̄ (true mass) ≈ 1.01–1.02; "steady to 2 %" holds only on the coarse
   grid and only to t ≈ 300.
8. **Lone stars are stationary** (phase5/spherical): an accurate 1+1 nonlinear evolution keeps both campaign
   stars (ε = ±1) stationary to t = 1500 (drifts ≤ 1e-6, 4th-order convergent; unstable-star growth rates match
   the radial analysis to < 0.6 %). Shirokov's lone-star drifts are 4–5 orders larger, so they are 3D-code or
   initial-data artefacts; his 5 %-radius seam is not the cause. The frequency→mass proxy holds to a few %.

Novelty status of each result: `notes/novelty/RESULTS.md`.

9. **No non-radial instability found** (phase2/nonradial): for the phantom star, l = 2, 3 polar perturbations have no
   growing mode (real γ up to 300 / 1000, covering the point-mass rates 137 / 686; no complex growing mode in
   0.1 ≤ Re s ≤ 2, |Im s| ≤ 2.5), and the axial sector is stable for all l. The interior supplies κ_eff > κ_c, so the
   point-mass instability is not triggered. Scans, not proofs: l ≥ 4, Re s < 0.1 and compact stars remain open.

## Setup
```
pip install sympy mpmath numpy scipy
python phase0/verify_zerilli.py   # several minutes (symbolic limits)
python phase0/mode_width.py
```
Each phase directory has a `RESULTS.md` and per-script `*_output.txt`.

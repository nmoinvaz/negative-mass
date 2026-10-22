# negative-mass

Working notes and code on negative (active) mass in general relativity: Bondi's 1957 runaway, the
stability of negative-mass Schwarzschild, regular negative-mass sources, and what remains open.

- `PLAN.md` — phases, status, next steps
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
   Exactly one unstable mode per l (checked l ≤ 8), localized at the singularity, growth ∝ l⁴.
2. **A regular star hides that instability.** φ₀ = 1/ψ₀ is an exact zero mode of the transformed operator, giving
   a closed-form threshold: with a non-energy-injecting surface condition (κ ≥ 0) a star is stable once
   R > r_*(l) (l = 2: 2.279|M|). The l⁴ growth belongs to the point singularity only.
3. **Shirokov's phantom boson stars reproduced.** Profiles match his tables to 1.4×10⁻⁵; his "ADM masses" are
   the mass enclosed where φ = 5 % φ_c and run 1.5–2.3 % below the true ADM mass.
4. **The phantom campaign star is radially stable.** The family loses stability exactly at dN/dω = 0; the
   turning-point criterion carries over to negative mass when stated in N (not M).
5. **The phantom point particle is Bondi's all-negative mass**: one parameter is inertial, passive and active
   mass; gravity alone cannot distinguish "positive inertial, negative active". At 1PN the dipole midpoint
   acceleration is exactly m/d² (no correction); a matched pair slowly contracts (d̈ = −E_tot/d²).
6. **Mass bounds:** per slice, M ≥ −(phantom's weighted energy), tight to 0.1 % for the phantom star; no bound
   uniform in time. Whether radiation can drive the Bondi mass to −∞ is open.
7. **Review of Shirokov's code:** the runaway physics is implemented correctly (signs, couplings, diagnostics
   verified). Two real bugs off the paper's path (boosted-star Π sign; Ψ4 matter term zeroed). Precision
   claims are weaker than stated: the a = M̄/d² agreement rests on the coarse grid (fine grid 7 % higher at
   d = 10), and the acceleration rises ~11 % over the run on both fine grids. That rise is not gauge, separation,
   constraint or boundary error and barely converges between the two fine grids; field-content growth runs
   the wrong way across grids. Only the stars' slowly drifting frequency tracks it (implying ~9–11 % mass growth),
   unconfirmed. Early-time fine-grid a·d²/M̄ (true mass) ≈ 1.01–1.02; "steady to 2 %" holds only on the coarse
   grid and only to t ≈ 300.

In progress: non-radial (l ≥ 2) stability of the phantom star.

## Setup
```
pip install sympy mpmath numpy scipy
python phase0/verify_zerilli.py   # several minutes (symbolic limits)
python phase0/mode_width.py
```
Each phase directory has a `RESULTS.md` and per-script `*_output.txt`.

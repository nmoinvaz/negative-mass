# Plan

## Phase 0 — Ground truth ✅
Reproduce GHI's q_C and check the published modes. Done: see `phase0/RESULTS.md`.
Outcome: negative-mass Schwarzschild is linearly unstable (Dotti–Gleiser 2008); GHI's stability rests on a
discontinuous metric perturbation. Gap 5 (Paranjape vs. Shirokov) also closed.

## Phase 1 — Stability map ✅ (superseded)
The boundary-condition ambiguity is resolved by the DG08 intertwiner. No further work unless we want an
independent reproduction of the DG08 spectrum (cheap with the Phase 0 code).

## Phase 2 — Regular negative-mass stars ✅ (see phase2/RESULTS.md; open: l ≥ 4, slow modes, compact stars)
- ✅ Phantom star radially stable; no l = 2, 3 polar or axial instability; κ_eff > κ_c (phase2/nonradial).
- ✅ Vacuum exterior cannot host an even-parity instability when |M|/R < 0.374 (Robin κ ≥ 0).
- ✅ Shirokov's canonical and phantom stars reproduced (profiles exact; paper's masses are enclosed at φ=5 % of φ_c, ours are true ADM).
- Build static, spherically symmetric negative-ADM-mass stars (Shirokov's sextic phantom soliton, and a
  perfect-fluid alternative) as an ODE shooting problem.
- Compute their even-parity linear spectrum for l = 0, 1, 2, 3, ….
- Test the heuristic: no unstable l ≥ 2 mode for k ≳ 8(|M|/R)²; find the critical compactness where one appears.
- Check convergence of the unstable spectrum to the Gleiser–Dotti modes as R → 0.
- Question: is the point-singularity problem ill-posed (growth unbounded in l), and does any regular star cure it?

## Phase 3 — Conservation and the third law ✅ (see phase3/RESULTS.md)
Result: the phantom particle is Bondi's all-negative mass (μ = s·m in every role; gravity alone cannot tell
them apart); Bergmann holds. 1PN: midpoint acceleration exactly m/d² (no 1PN correction); proper separation
obeys d̈ = −E_tot/d², rigid only if |μ₋| = μ₊(1+μ₊/d). Conserved E, signed P, J, G verified. 1PN good to ~1%
up to v≈0.4. Shirokov's +6–8% excess is not a 1PN or finite-size effect.

(original plan:)
First post-Newtonian two-body expansion of G = 8π(T₊ − T₋): which mass combinations appear, and does
Bergmann's "active ≠ passive is impossible" survive?

## Phase 4 — Bondi mass bound ✅ (see phase4/RESULTS.md)
Result: per-slice bound M ≥ −W₋ (Witten; exact identity for conformally flat data) and an a-priori
M ≥ −E₋/(1−κ) for κ<1, tight to 0.1% for the phantom star; no uniform-in-time bound (E₋, W₋ not conserved;
no universal constant). GW luminosity computed (phase4/radiation): L₀ = (24/35)m⁶/d⁶ in the rest frame; M_B → −∞ not established
(toy back-reaction model suggests a finite loss). No prior bound
of this kind found (closest: Shi–Tam 2004).

(original plan:)
Witten spinor argument with a signed source: lower bound on the Bondi mass, or a construction showing
unbounded energy extraction. Stop after ~2 weeks if it trivializes.

## Phase 5 — Nonlinear checks
✅ 1+1 lone-star evolution (phase5/spherical/RESULTS.md): both campaign stars stationary to t=1500; Shirokov's
lone-star drifts not physical. Remaining: 3D items below (need GPUs).

1+1 evolution of the lone phantom star (is the 0.63%/100 decay physical?); 3+1 Shirokov N = 192 test and
an off-axis kick.

## Later — Ghost-wormhole threshold vs. Choptuik critical collapse (to investigate)
From the literature check in `notes/critical_ghost.md`: the ghost-wormhole threshold is type I (finite black-hole
mass, measured by Shinkai–Hayward 2002) and Choptuik's is type II. No paper connects the two, and in the normal + ghost
model they sit in disjoint sets of initial data. Saved to investigate later:
- **Where a real connection could still exist** (none studied):
  1. self-similar solutions in which the normal and ghost fields boost into each other (the Lorentzian analogue of
     Hirschmann–Eardley's rotating complex-scalar solution);
  2. Brans–Dicke gravity with ω < −3/2 (a ghost scalar in the Einstein frame);
  3. the outcome map for regular-centre normal + ghost initial data, where temporary throats have been seen.
- **Untested prediction:** a small ghost admixture leaves Choptuik's critical solution unchanged (same Δ and γ).
  Reason: its ghost-sector perturbations are the same problem as a complex scalar's imaginary-part perturbations, which
  all decay except one symmetry mode. Test: compute those perturbations of the choptuik-large-d solution across D and
  look for an eigenvalue crossing zero (that would signal a genuinely two-field branch).
- **Can our solvers represent a wormhole throat?** Not as they stand: choptuik-large-d's areal radial coordinate (and
  this repo's 1+1 code) cannot pass a throat, which shows up like an apparent horizon. Investigate whether a change of
  coordinates (proper radius, or double-null/characteristic coordinates as Shinkai–Hayward used) makes it possible,
  and what it would cost.

## Out of scope
Quantum vacuum instability of the phantom sector.

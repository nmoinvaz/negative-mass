# Plan

## Phase 0 — Ground truth ✅
Reproduce GHI's q_C and check the published modes. Done: see `phase0/RESULTS.md`.
Outcome: negative-mass Schwarzschild is linearly unstable (Dotti–Gleiser 2008); GHI's stability rests on a
discontinuous metric perturbation. Gap 5 (Paranjape vs. Shirokov) also closed.

## Phase 1 — Stability map ✅ (superseded)
The boundary-condition ambiguity is resolved by the DG08 intertwiner. No further work unless we want an
independent reproduction of the DG08 spectrum (cheap with the Phase 0 code).

## Phase 2 — Regular negative-mass stars (in progress — see phase2/RESULTS.md)
- ✅ Vacuum exterior cannot host an even-parity instability when |M|/R < 0.374 (Robin κ ≥ 0).
- ✅ Shirokov's canonical and phantom stars reproduced (radius exact; masses +2 % common factor, open).
- Build static, spherically symmetric negative-ADM-mass stars (Shirokov's sextic phantom soliton, and a
  perfect-fluid alternative) as an ODE shooting problem.
- Compute their even-parity linear spectrum for l = 0, 1, 2, 3, ….
- Test the heuristic: no unstable l ≥ 2 mode for k ≳ 8(|M|/R)²; find the critical compactness where one appears.
- Check convergence of the unstable spectrum to the Gleiser–Dotti modes as R → 0.
- Question: is the point-singularity problem ill-posed (growth unbounded in l), and does any regular star cure it?

## Phase 3 — Conservation and the third law
First post-Newtonian two-body expansion of G = 8π(T₊ − T₋): which mass combinations appear, and does
Bergmann's "active ≠ passive is impossible" survive?

## Phase 4 — Bondi mass bound (high risk)
Witten spinor argument with a signed source: lower bound on the Bondi mass, or a construction showing
unbounded energy extraction. Stop after ~2 weeks if it trivializes.

## Phase 5 — Nonlinear checks (optional)
1+1 evolution of the lone phantom star (is the 0.63%/100 decay physical?); 3+1 Shirokov N = 192 test and
an off-axis kick.

## Out of scope
Quantum vacuum instability of the phantom sector.

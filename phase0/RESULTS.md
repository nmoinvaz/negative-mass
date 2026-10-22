# Phase 0 results — negative-mass Schwarzschild, even-parity perturbations

Scripts: `verify_zerilli.py` (output in `verify_zerilli_output.txt`), `mode_width.py`.
Conventions follow Gleiser–Dotti: f = 1 − 2M/r, dr*/dr = 1/f, λ = (l−1)(l+2)/2, M < 0, μ = −2M.

## Verified (symbolic residuals exactly 0)

- Both pieces of the GHI static l = 2 solution (their eq. 19) solve the Zerilli equation.
- The GHI combination C3 = C4[log(7/3) − 4/9] (μ = 1) has no pole at r_C = 3μ/4 and vanishes there as (r − r_C)².
- Fitting Φ ≈ a₀ r*^{1/2} log r* + b₀ r*^{1/2} near r = 0 gives q = b₀/a₀ → 4.3877427 (r ~ 1e-10),
  vs. GHI's q_C = 2 + log(98/9) = 4.3877429.
- The Gleiser–Dotti unstable mode Ψ = r(r−2M)^k e^{k(t−r)/2|M|} / (2λr + 6M) solves Zerilli for l = 2…5,
  growth rate k/(2|M|), k = (l−1)l(l+1)(l+2)/6.
- DG08's intertwiner seed ψ₀ (l = 2) solves the static Zerilli equation.

## The GHI stability argument fails at r_C, not at the singularity

GHI truncate their marginal mode (Φ = 0 for r > r_C). Reconstructing the metric functions from it
(μ = 1, l = 2): K → 0 but dK/dr → 4.063, and **H₂ → 3.048 from the left vs. 0 from the right**.
The implied metric perturbation is discontinuous, confirming Gleiser–Dotti's objection.

Smooth metric perturbations instead have a simple pole in Ψ at r_s = −3M/λ with
c₀ = λ² c₋₁ / (3M(3 + 2λ)) (DG08 Lemma 1). The unstable mode satisfies it:
c₀ = −0.533655 computed vs. −0.533655 predicted (M = −1, l = 2).

**Status of gap 4:** closed in the literature. DG08 remove the pole with an intertwining operator; the
transformed operator has a unique self-adjoint extension matching the physically motivated r = 0
condition, and generic data supported away from the singularity excites the unstable mode.
Negative-mass Schwarzschild is linearly unstable.

## New observation: the instability lives at the singularity

| l | growth × \|M\| | 1/e radius / \|M\| | 2√(2/k) |
|---|---|---|---|
| 2 | 2 | 1.766 | 1.414 |
| 3 | 10 | 0.701 | 0.632 |
| 4 | 30 | 0.388 | 0.365 |
| 5 | 70 | 0.249 | 0.239 |
| 6 | 140 | 0.174 | 0.169 |
| 7 | 252 | 0.129 | 0.126 |

Width shrinks ~ l⁻², growth rises ~ l⁴. **Heuristic (untested):** a regular star of radius R removes the
modes with width ≲ R, i.e. all l ≥ 2 once k ≳ 8(|M|/R)². For weak compactness (Shirokov's star,
|M|/R ≈ 3×10⁻³) no singular mode survives outside the matter — consistent with his lone phantom
star lasting to t = 1000. Open: the growth rate is unbounded in l, which for the point singularity
looks closer to ill-posedness than an ordinary instability.

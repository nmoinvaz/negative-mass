# The ghost-field wormhole's instability at large D

October 2026. A follow-up to [`../critical_ghost.md`](../critical_ghost.md) (CG) §1.1, which noted from Torii and
Shinkai's table [TS13] that the wormhole's single unstable mode grows faster in higher dimensions, roughly like
√(2D)/a. This note derives that dependence.

**The system.** The static, massless, symmetric wormhole held open by a ghost scalar (Ellis–Bronnikov) in D spacetime
dimensions, with throat radius a:

```math
ds^2 = -dt^2 + dl^2 + r(l)^2\,d\Omega_{D-2}^2, \qquad r'^2 = 1 - (a/r)^{2(D-3)} .
```

Its radial perturbations have exactly one growing mode, e^{Wt} [GGS09a for D = 4; TS13 numerically for D = 4–11, 20,
50, 100]. Since g_tt = −1, t is proper time everywhere, and 1/W is the e-folding time of the throat's collapse or
expansion.

Status words: **derived** (closed-form matched asymptotics, each step checked numerically), **computed** (numerics,
converged), **identified** (a closed form found by an integer-relation search on a computed number).

## Result

```math
W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\pi^2/3 + 4\ln^2 2}{D-3} + O(D^{-2}) .
```

- The first two terms are **derived**.
- The third is **computed** to ten digits from the third-order matching integral, Ω₁ = 10.4233603788. It is
  **identified** as 2π²/3 + 8 ln²2 = 10.42336037874, with 3Ω₁ − 2π² − 24 ln²2 = 0 found by PSLQ at tolerance 10⁻⁹.
- Against the full numerics the three terms leave about 73/(2(D − 3))² at D = 200 and 500, as an O(1/D²) remainder
  should. D = 1000 gives 66 ± 6, within that value's numerical error (table below).

**What it means.**
- **Faster in higher D.** The e-folding time is a/√(2(D − 3)) at large D. The exact values are 0.846 a (D = 4),
  0.579 a (D = 5) and 0.462 a (D = 6), falling to 0.072 a at D = 100.
- **Where the mode lives (derived).** Near the throat, in the shell of width about a/D where the geometry differs from
  flat space, the throat moves rigidly: δr is constant there, with δr′/δr = −1/a at its outer edge. Outside the shell
  space is nearly flat, and the mode is a Gaussian, δr ∝ (r/a) e^{−r²/a²}, to leading order. So at large D the
  instability is a flat-space mode with the throat acting as a boundary.
- **A prediction for threshold evolutions.** CG §1.3 reads the knife edge between collapse and expansion as type I
  critical behaviour, with the static wormhole as critical solution. Then near the threshold the time before the
  wormhole collapses or expands grows like T ≈ (1/W) ln(1/|p − p*|). The slope 1/W should be 0.579 a at D = 5 and 0.462 a
  at D = 6. Shinkai and Torii evolved D = 5 and 6 [ST17] but did not measure this slope. At D = 4 it is the 0.846 a that
  Shinkai and Hayward's 2002 runs already show (CG §1.3).

## Method

### 1. The master equation

Gauge δφ = 0, which is regular everywhere because φ′ never vanishes. The perturbations are ν = N e^{Wt}, λ = L e^{Wt}
and r = r₀ + Y e^{Wt} in ds² = −e^{2ν}dt² + e^{2λ}dl² + r²dΩ². With n = D − 2 and a = 1 (`derive.py`):
- the tl equation gives L = Y′/r₀′;
- the scalar equation gives N = L − nY/r₀ (an integration function of t is removed by reparametrizing t);
- the angular equation gives

```math
Y'' - \left(\frac{2r_0''}{r_0'} + \frac{n\,r_0'}{r_0}\right)Y' + \left(\frac{n\,r_0'^2}{r_0^2} + \frac{r_0''}{r_0}\right)Y = W^2Y .
```

The tt and ll equations then hold identically (`master.py check`, symbolic). Its control, a wrong coefficient of Y′,
leaves nonzero residuals. The throat (r₀′ = 0) is a regular singular point with exponents 0 and 3. The mode that moves
the throat is the even one, with Y(0) = 1.

### 2. Numerics

Two independent integrations:
- in proper distance l (`master.py modes`, for D ≤ 10);
- in the areal radius r (`modes_r.py`, for D ≥ 20). It uses a Frobenius start at the throat and matches to the exact
  far-field solution r^{(D−1)/2} K_{(D−3)/2}(W r).

Both find exactly one growing mode, and they agree with [TS13]'s Table I. TS13 used a different master variable, a
gauge-invariant combination regularized at the throat.

| D | W²a² here | [TS13] | 2(D − 3) − 4 ln 2 + Ω₁/(2(D − 3)) |
|---|---|---|---|
| 4 | 1.39708 | 1.39705243 | (not asymptotic) |
| 10 | 12.0455 | 12.0442650 | 11.972 |
| 20 | 31.5751101276 | 31.5751101285 | 31.53398 |
| 50 | 91.3457759034 | 91.3457759137 | 91.33830 |
| 100 | 191.2830176842 | 191.2830177297 | 191.28114 |
| 200 | 391.2543342170 | | 391.25387 |
| 500 | 991.2379714266 | | 991.23790 |
| 1000 | 1991.2326553188 | | 1991.23264 |

- The l-integration is accurate to about 10⁻³ at D ≤ 10. At large D its start near the throat drifts, which is why the
  r-integration exists.
- The r-integration's values move by 4×10⁻¹⁰ to 1.5×10⁻⁶ from D = 20 to 1000 when its numerical parameters are changed
  (`logs_modes_r.txt`).
- A fit of the r-integration values in 1/D gives W²a² − 2(D − 2) = −4.77259 ± 10⁻⁵, against −2 − 4 ln 2 = −4.7725887.

### 3. Large D: matched asymptotics

**Variables.** Near the throat, R = (r/a)^s with s = 2(D − 3), ε = 1/s, ℓ = ln R and r = e^{εℓ}. Write
W² = s + Ω with Ω = Ω₀ + εΩ₁ + …. The master equation becomes exactly

```math
Y_R = \varepsilon\,\frac{L\,r}{R}, \qquad (R-1)L_R = (1+\varepsilon\Omega)\,rY + \frac{L}{2R} + \left(\tfrac12+\varepsilon\right)\frac{L(R-1)}{R} - \left(\tfrac12+\varepsilon\right)\frac{Y(R-1)}{Rr} - \frac{Y}{2rR},
```

where L = dY/dr. At each order in ε this gives (R − 1)L_k′ − L_k/2 = h_k(R) (`matched.py`).

**The condition at each order.** The solution analytic at the throat is

  L_k = −2h_k(1) + (R − 1)^{1/2} ∫₁^R (h_k − h_k(1))(R′ − 1)^{−3/2} dR′.

At large R it contains G_k R^{1/2}, with G_k = ∫₁^∞ (h_k − h_k(1))(R − 1)^{−3/2} dR. But R^{1/2} = (r/a)^{D−3} is the
far zone's *growing* solution. So the mode requires G_k = 0 at every order.

**The orders** (closed forms; u = √(R − 1)):
- **Order 0:** Y₀ = 1, L₀ = −1.
- **Order 1:**
  - h₁ = Ω₀ + ln R − 2(R − 1)/R.
  - G₁ = 0 identically (∫ ln R (R − 1)^{−3/2} = 2π = ∫ 2(R − 1)/R · (R − 1)^{−3/2}).
  - L₁ = −2Ω₀ − 2 ln R exactly, and Y₁ = −ln R.
- **Order 2:**
  - Y₂ = −2Ω₀ ln R − (3/2)ln²R.
  - h₂ = Ω₁ − 2Ω₀(R − 1)/R − Ω₀ ln R − 2 ln²R.
  - G₂ = −4πΩ₀ − 2J, with J = ∫₁^∞ ln²R (R − 1)^{−3/2} dR = 8∫₀^∞ ln(1 + u²)/(1 + u²) du = 8π ln 2.
  - So **Ω₀ = −4 ln 2.**
- **Order 3:**
  - L₂ = −2Ω₁ + 2Ω₀ ln R + 4 ln²R − 8Ω₀ u arctan u − 16u K(u), with K(u) = ∫₀^u ln(1 + t²)/(1 + t²) dt. Its R^{1/2}
    term vanishes exactly when Ω₀ = −4 ln 2, as it must.
  - G₃(Ω₁) is linear, with slope −4π (computed: −12.5663706143, against −4π = −12.5663706144). G₃ = 0 gives
    **Ω₁ = 10.4233603788** (`order3.py`; the two finer of three grids agree to 10⁻⁹).

**Matching.** The decaying far solution r^{(n+1)/2}K_ν(Wr), expanded in the overlap 1 ≪ R ≪ e^s, gives Y′/Y =
−1 − (2Ω₀ + 3 ln R)/s + …. The near zone gives the same, so the O(1/s) matching holds automatically. The constraint
comes only from the absence of the growing term.

## Limits

- Radial (spherically symmetric) perturbations only.
- The massless symmetric member only. The asymmetric (two masses), charged and massive-potential families are not done.
- The third term's closed form is an identification, not a derivation.
- **Prior work** (searched only briefly, three web searches). [TS13] give the numbers and say the time scale "becomes
  shorter as n becomes large", but give no formula. No large-D analysis of these wormholes was found. This is not a
  full literature check, and nothing here claims novelty.

## Files

- `derive.py`: the first-order equations.
- `master.py`: the symbolic check and its control (`logs_check.txt`), and the l-integration (`logs_modes_lowD.txt`).
- `modes_r.py`: the r-integration with its convergence check (`logs_modes_r.txt`).
- `matched.py`: the h_k from the series (`logs_matched.txt`).
- `order3.py`: Ω₁ (`logs_order3.txt`).

## References

- [TS13] T. Torii, H. Shinkai, "Wormholes in higher dimensional space-time: exact solutions and their linear stability
  analysis", PRD 88, 064027 (2013), arXiv:1309.2058. Table I read for this note (full-precision values quoted above).
- [GGS09a], [SH02], [ST17]: see CG's reference list.

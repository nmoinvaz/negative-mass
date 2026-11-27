# Polar modes of the phantom black universes and radial modes of the M-AdS wormholes

November 2026, #218 (branch `wip/i218-polar-mads`; raw `I218-dh610` in negative-mass-raw). Code `polar_derive.py` (derivation and its checks,
`polar_derive_output.txt`) and `polar.py` (numerics, `polar_output.txt`). Extends #212 (`README.md` in this directory) and `../critical_ghost.md` §1.1.
Status words as in AGENTS.md: counts and rates are *numerical evidence*; the equations are *checked symbolically* (sympy, in `polar_derive.py`); Sec. 1
holds one analytic argument at the Ellis wormhole (a positive potential), also *checked symbolically*.

## Result

- **No growing polar mode.** At the 20 black universes of #212 (c = −10 … −0.01, including the isolated member c = −1) and at the Ellis wormhole, polar ℓ = 1 and ℓ = 2 have no growing mode with real
  W in [0.01, 8] per b: 84 runs, R1 and R2 agree (the determinants at five W to 5.1e-9), no sign change and no interior minimum of the determinant. At the Ellis wormhole the absence is also proved,
  because the master potentials are positive: `L/R² − 3/R⁴` for ℓ ≥ 2 and `1/R⁴` for ℓ = 1 (Sec. 1). *Numerical evidence* for the members, *checked symbolically* at c = 0. Complex W is not searched (Sec. 7).
- **M-AdS wormholes: exactly one growing radial mode** at each of 12 values of c = 1e-4 … 10 (Dirichlet at the AdS end), W = 1.1820 … 15.221 per b, e-folding time 0.846 … 0.0657 b; R1 and R2 agree to 1.6e-10;
  `W²(c) → 1.39705243433` (the Ellis value) as c → 0⁺ (extrapolated intercept 3.1e-10), and `W²(c) − W²(0) ≈ 1.863 c + 1.30 c²`. This agrees in kind with [BKZ12] ("unstable", by regularized
  potentials and time-domain runs) and adds the count and the rates. *Numerical evidence.* A Chebyshev eigenvalue solver sharing none of the AdS series, throat crossing or matching gives the same counts and the same
  rates (S-B: 6.8e-5 at c = 0.01, 1.7e-6 or better for 0.03 ≤ c ≤ 5); one declared tolerance is **not met**: at c = 10 the finer collocation (N = 360, X = 60) differs by 5.9e-4 against the declared 2e-4
  (the coarser one, N = 240, X = 40, differs by 2.9e-4 against 2e-3).
- **Planted cases pass.** The derived polar system at c = 0 gives ℓ = 2, 3, 4 quasinormal modes 1.737692 − 0.305138i, 2.952415 − 0.409984i, 4.076253 − 0.449105i against [ABCKY23]'s 1.74 − 0.305i, 2.95 − 0.410i,
  4.08 − 0.449i; its ℓ = 1 mode, 0.681367187 − 0.617753448i, equals the printed master equation's at L = 2 (1.5e-12; 2.8e-11 against an eigenvalue code).
- **The printed text disagrees with the derivation in three places** (Sec. 5): (47), (69) and the statement that the ℓ = 1 polar mode is the scalar branch of [BCK18]. The ℓ = 1 problem is not published as a number; the ruling of the coordinator
  fixed the planted case (Sec. 3).
- 11 + 9 `--mutate=` controls, each computed in a development run before it was declared, all fire (Sec. 4).

## 1. The linearization

The family is that of #212 (b = 1, R² = ρ² + 1, φ₀′ = √2/R², ε = −1): `ds² = −A dt² + dρ²/A + R² dΩ²`, with `A″ = 2(A − 1)/R²` (two constants, c and m; 3m = −2c/π
for asymptotic flatness at ρ → +∞) and, from the tt and scalar equations,

```math
V = -\frac{A - 1 + \rho A'}{R^2},\qquad V_\varphi = -\frac{\sqrt2\,A'}{R^2},\qquad V_{\varphi\varphi} = -A'' + \frac{2\rho A'}{R^2}.
```

`R_ab − ε ∂_aφ ∂_bφ − V g_ab = 0` and `ε □φ − V_φ = 0` are linearized about it with perturbations `e^{Wt} Y_l(θ)` (`polar_derive.py`: `linearize`). The two checks that decide
the linearization itself are D1 (background residual 0; V against `stability.family_mp` to 1.5e-40) and D2 (the Lie derivative of g and φ₀ along a generic
vector ξ = (T Y, X Y, S ∂_θY) satisfies all eight linearized equations identically; with `lin_sign` mutated it does not).

**l ≥ 2, Regge–Wheeler gauge** (`h_tA = h_ρA = G = 0`): `h_tt = A H0 Y`, `h_tρ = H1 Y`, `h_ρρ = H2 Y/A`, `h_θθ = R² K Y`, `δφ = F Y`. This gauge is complete for l ≥ 2 and its
equations contain only `1/A` and `ρ/R²`, never `1/R′`: no singular point at the throat. Decisive steps (D3): the traceless angular equation is `H2 = H0`; the tθ, tρ, ρθ equations give
`H1′, K′, H0′`; the tt, ρρ and trace-angular equations, after those, are one algebraic relation (ratios `1/A²` and `−2R²/A`), which gives `F′`; `d/dρ` of that `F′` equals `F″` of the scalar
equation and the four second-order equations vanish on the result. With `Π = H1/W`, `Φ = (H0, Π, K, F)`:

```math
\Phi' = (M_0 + W^2 M_1)\,\Phi,\qquad \text{rows } (H_0', \Pi', K', F'),\ R^2=\rho^2+1,
```
```
H0': [ (−A₁R² + Aρ)/(A R²),            (A L/2 + W²R²)/(A R²),            (A₁R²/2 − Aρ)/(A R²),                    −√2/R² ]
Π' : [ 1/A,                             −A₁/A,                            1/A,                                      0 ]
K' : [ ρ/R²,                            L/(2R²),                          (A₁R²/2 − Aρ)/(A R²),                     √2/R² ]
F' : [ √2(−A₁ρ + 2A − L)/(4A),          √2(−A₁L − 4W²ρ)/(8A),             √2(−A₁²R² + 2A₁Aρ + 2AL − 4A + 4W²R²)/(8A²),   (A₁R² − 4Aρ)/(2AR²) ]
```

with `A₁ = A′`, `L = l(l+1)`. W enters only through W²: the system is time-reversal symmetric, and the eigenvalue problem is the linear pencil above. (`polar.py` carries these as
`m4_entries`, emitted by `polar_derive.py emit`; D7 checks the table against the derivation.)

**l = 1.** There is no traceless equation and the RW gauge leaves one function of ρ free. The gauge `h_tA = h_ρA = 0, h_tρ = 0` (the gauge `H1 = 0` of [ABCKY23] Sec. 3.3, whose `K = −H2` it shares) gives `K = −H2` and a
3-system for `(H0, H2, F)`:

```
H0': [ (−A₁R² + 2Aρ)/(2AR²),   (−A₁R² + Aρ)/(AR²),   −√2/R² ]
H2': [ 0,                       (A₁R² − 4Aρ)/(2AR²),   −√2/R² ]
F' : [ √2(−A₁ρR² + 2Aρ² − 2R²)/(4AR²),   √2(A₁²(ρ⁴+1)/8 − A₁ρ(−A₁ρ + AR²)/4 + A²/2 − W²R⁴/2)/(A²R²),   (A₁R² − 4Aρ)/(2AR²) ]
```

(D6a, b: consistency and the second-order equations.) The one residual gauge freedom is `ξ^θ = c(ρ)` with `c = √A/R`, `ξ^ρ = −A R² c′`, `ξ^t = R² W c/A`; its fields `g = (H0, H2, F)` are in closed form
(`gauge_l1`) and solve the 3-system exactly (D6c, residual 0). A solution that decays at both ends modulo this gauge solution is `det[L | R | g] = 0`, L and R the solutions that are ingoing at the left end and decaying
at the right end.

**Ellis wormhole (A = 1).** (D4a) the l ≥ 2 system reproduces [ABCKY23] (44), (45), (46) exactly (u = F/√2, ω = iW); (D5a, b) F decouples, `F″ = −(2ρ/R²) F′ + b F`, and Z = R F obeys

```math
Z'' = (V_l + W^2)\,Z,\qquad V_l = \frac{L}{R^2} - \frac{3}{R^4},
```

their (52), (53): the scalar and the gravitational polar fields have the same master equation and the polar spectrum is two-fold degenerate. *Consequence* (checked symbolically): for L ≥ 6, `V_l > 0`
(`L R² − 3 > 0`), so `−Z″ + V Z = −W² Z` has no solution with W² > 0: **no growing polar l ≥ 2 mode at the Ellis wormhole**, for every l.
(D6d) for l = 1 the gauge-invariant scalar `χ = F − (g_F/g_H2) H2` obeys `χ″ + (2ρ/R²)χ′ − W²χ = 0`, i.e. `Z = Rχ`: `Z″ = (1/R⁴ + W²)Z`, `V = 1/R⁴`. **No growing l = 1 mode at the Ellis
wormhole** either (V > 0). (D6e) With `w = −ρ/R²`: `w² − w′ = 1/R⁴` and `w² + w′ = 2/R² − 3/R⁴` = `V_l` at L = 2: the two are Darboux partners, so the l = 1 problem is isospectral to the
printed l ≥ 2 master equation extended to l = 1; (D6f) w is regular on the real line (poles at ρ = ±i) and tends to 0 at both ends, so `Q = d/dρ + w` maps `e^{∓Wρ}` to `(∓W + w)e^{∓Wρ}` and the
quasinormal-mode conditions are kept.

## 2. Boundary conditions

- **Right end** (all members): both fields decay, `e^{−Wρ}`; the frame is integrated inward from `ρ_far = max(rmax0, rmaxW/W)`, where the decaying solutions dominate, with Gram–Schmidt (QR) between chunks.
- **Horizon** (c < 0, ρ_h < 0 for c > −1, ρ_h = 0 at c = −1, ρ_h > 0 below): ingoing. The RW fields are singular at a horizon, so the system is made Fuchsian first:
  `F = F̃/z` (l = 2), `H2 = z H̃2` (l = 1), `z = ρ − ρ_h`. The residue matrix then has the exponents
  `{1 + μ, −1 + μ}` (ingoing) and `{1 − μ, −1 − μ}` (outgoing), μ = W/A′(ρ_h) (l = 2), and `{μ}` (ingoing), `{−μ, −1/2}` (outgoing, gauge) for l = 1; ingoing is the family
  whose exponents rise with W, the regular solution of the ingoing Eddington–Finkelstein form (`e^{Wv}`, `v = t + ∫dρ/A`). The ingoing solutions are Frobenius series, with Taylor coefficients of
  `z M′(z)` taken by FFT on a circle of radius 0.35; the recursion at the resonance `s₂ + 2 = s₁` is consistent (no logarithm; the worst relative inconsistency is checked < 1e-6 in every run).
  At small W the second ingoing exponent `−1 + μ` is below the outgoing `1 − μ`, so a bare leading-order start would be contaminated; the series (10 and 12 terms) removes that.
- **Ellis** (c = 0): the second infinity, decaying `e^{+Wρ}`, ρ → −∞.
- **M-AdS wormholes** (c > 0, radial): A ~ 2cρ² as ρ → −∞, `Y → α + βρ` in the gauge δφ = 0; Dirichlet (no source, `Ψ = −√2 Y/ρ → 0`) is β = 0, the series
  `Y = 1 + y₂/ρ² + y₄/ρ⁴ + …`, `y₂ = (W² − 4c² − 4c)/(24c²)`, `y₄ = (W⁴ − 40W²c + 272c⁴ + 416c³ + 24c²(4 − 3W²))/(1920c⁴)` (sympy), started at `|ρ| ≥ max(100, 400 W/c, 100/√c)`,
  continued to ρ = −δ, through the throat by the complex arc of #212, and matched at ρ = 1 to the solution of #212 that decays at ρ → +∞.

## 3. Method, checks and planted cases

Shooting at real W, as in #212 but on a 2-frame (l = 2) or 1-frame plus gauge column (l = 1) determinant; the determinant at the matching point is `D(W)`; a growing mode is a root. Real roots are sign changes of `D` on a
log grid in W (60 points at R1, 150 at R2; Brent), and, because the Ellis polar fields are degenerate (two identical master equations) and its roots are double, interior minima of `|D|` below 1e-3 of both neighbors,
refined by a bounded search and counted twice. The range searched is W in [0.01, 8] per b for the polar modes and [0.01, 24] for the M-AdS modes. R1: rtol 1e-9, far boundary max(150, 30/W);
R2: rtol 1e-12, max(300, 45/W) (the full settings are in the raw README).

| check | what | result |
|---|---|---|
| D1, D1b | background on the family, sympy; V against `stability.family_mp` | residual 0; 1.5e-40 |
| D2 | pure-gauge perturbation satisfies the eight linearized equations | residuals 0 |
| D3a–e | RW reduction, l ≥ 2 | consistent; matrix linear in W² |
| D4a, b | l ≥ 2 system against [ABCKY23] (44)–(47) | (44)–(46) equal; (47) differs by K(ℓ+2)(ℓ−1)/2 |
| D5a, b | F decouples at c = 0; master potential | `L/R² − 3/R⁴`, [ABCKY23] (53) |
| D6a–f | l = 1: consistency, gauge solution, χ-equation, Darboux relation with w regular and → 0 | residuals 0; `V = 1/R⁴`; `w² ∓ w′` |
| D7 | coefficient table in `polar.py` against the derivation | 2.8e-16 |
| P2 | Ellis polar QNM, ℓ = 2, 3, 4, against [ABCKY23] Table 1 (M = 0 row) | 1.737692 − 0.305138i, 2.952415 − 0.409984i, 4.076253 − 0.449105i against 1.74 − 0.305i, 2.95 − 0.410i, 4.08 − 0.449i |
| P3a, b | Ellis polar ℓ = 1 QNM against the master equation at L = 2 (shooting; complex-scaled eigenvalue problem) | 0.681367187 − 0.617753448i; 1.5e-12, 2.8e-11 |
| P4 | Ellis at L = −1: double roots of D against finite differences | 0.421976, 1.473675 |
| P5 | M-AdS radial mode as c → 0⁺ | `W²(c) − 1.39705243433350684 = 1.8628634e-4, 1.864029e-3, 1.875713e-2` at c = 1e-4, 1e-3, 1e-2; intercept 3.1e-10 |

No number is published for the polar ℓ = 1 spectrum of the Ellis wormhole: [ABCKY23] give the equation (69) and refer to [BCK18], which has figures only. By the coordinator's ruling (hc, 2026-11-27) the planted case for ℓ = 1 is P3a and P3b
with the symbolic checks D6c–f (the pure-gauge solution; the Darboux relation, with w regular and → 0) and the control `l1_pot` (the test-scalar potential 2/R² + 1/R⁴ in its place, which misses P3a and P3b). ℓ = 2 keeps [ABCKY23] Table 1.
The determinant has double roots at the Ellis wormhole (P4 plants them), so the count of record there is also checked by the positivity argument of Sec. 1.

## 4. The counts

W is in units 1/b. "Modes" is `n_roots`: sign changes plus twice the interior minima of |D| below 1e-3 of both neighbors. The last column is the largest relative difference of the determinant between R1 and R2
at W = 0.01, 0.1, 1, 3, 8; the criterion was 1e-6. `polar_output.txt` has the same tables with the status column.

| c | ρ_h | κ | modes ℓ = 2 (R1, R2) | modes ℓ = 1 (R1, R2) | min \|D\| ℓ = 2 | min \|D\| ℓ = 1 | max rel. diff of D, R1 v R2 (ℓ = 2, ℓ = 1) |
|---|---|---|---|---|---|---|---|
| -10 | 4.197082 | 0.116522 | 0, 0 | 0, 0 | 5.88e-04 | 2.39e-03 | 6.4e-11, 1.4e-09 |
| -3 | 1.118725 | 0.351374 | 0, 0 | 0, 0 | 6.71e-03 | 2.73e-02 | 7.3e-12, 8.7e-10 |
| -2 | 0.621184 | 0.470502 | 0, 0 | 0, 0 | 1.22e-02 | 2.82e-02 | 9.9e-12, 5.3e-10 |
| -1.5 | 0.339617 | 0.551682 | 0, 0 | 0, 0 | 1.97e-02 | 1.62e-02 | 1.3e-11, 1.3e-09 |
| -1.1 | 0.075763 | 0.620954 | 0, 0 | 0, 0 | 2.60e-02 | 1.71e-02 | 1.5e-11, 7.8e-10 |
| -1.01 | 0.007824 | 0.635123 | 0, 0 | 0, 0 | 2.74e-02 | 1.72e-02 | 1.1e-11, 2.9e-10 |
| -1.001 | 0.000785 | 0.636471 | 0, 0 | 0, 0 | 2.75e-02 | 1.72e-02 | 1.1e-11, 3.3e-10 |
| -1.0001 | 0.000079 | 0.636605 | 0, 0 | 0, 0 | 2.75e-02 | 1.72e-02 | 1.1e-11, 3.3e-10 |
| -1 | 0.000000 | 0.636620 | 0, 0 | 0, 0 | 2.75e-02 | 1.72e-02 | 1.1e-11, 3.3e-10 |
| -0.9999 | -0.000079 | 0.636635 | 0, 0 | 0, 0 | 2.75e-02 | 1.72e-02 | 1.1e-11, 3.3e-10 |
| -0.999 | -0.000786 | 0.636768 | 0, 0 | 0, 0 | 2.75e-02 | 1.72e-02 | 1.1e-11, 3.4e-10 |
| -0.99 | -0.007884 | 0.638098 | 0, 0 | 0, 0 | 2.77e-02 | 1.72e-02 | 1.1e-11, 3.6e-10 |
| -0.95 | -0.040055 | 0.643811 | 0, 0 | 0, 0 | 2.83e-02 | 1.72e-02 | 1.1e-11, 4.7e-10 |
| -0.9 | -0.081826 | 0.650428 | 0, 0 | 0, 0 | 2.91e-02 | 1.71e-02 | 1.2e-11, 6.7e-10 |
| -0.7 | -0.271920 | 0.668151 | 0, 0 | 0, 0 | 3.22e-02 | 1.65e-02 | 9.8e-12, 1.9e-10 |
| -0.5 | -0.524073 | 0.660873 | 0, 0 | 0, 0 | 3.53e-02 | 1.50e-02 | 1.4e-11, 4.2e-10 |
| -0.3 | -0.926754 | 0.601300 | 0, 0 | 0, 0 | 3.84e-02 | 1.18e-02 | 1.8e-11, 5.5e-10 |
| -0.2 | -1.284188 | 0.532819 | 0, 0 | 0, 0 | 3.99e-02 | 9.27e-03 | 1.8e-11, 5.2e-10 |
| -0.1 | -2.024885 | 0.409506 | 0, 0 | 0, 0 | 4.14e-02 | 5.78e-03 | 1.6e-11, 6.5e-10 |
| -0.01 | -7.002156 | 0.140086 | 0, 0 | 0, 0 | 4.28e-02 | 1.52e-03 | 7.2e-12, 7.9e-10 |
| 0 | - | - | 0, 0 | 0, 0 | 8.30e-02 | 2.19e-04 | 1.1e-11, 5.1e-09 |

M-AdS wormholes, W in [0.01, 24] per b, Dirichlet at the AdS end; the criteria were 1e-6 (R1 v R2) and 1e-6 (monodromy, the relative imaginary part of Y after the throat arc, so no logarithm):

| c | modes (R1, R2) | W per b (R2) | W² | e-folding time 1/W | R1 v R2 | monodromy |
|---|---|---|---|---|---|---|
| 0.0001 | 1, 1 | 1.1820485273 | 1.39723872 | 0.84599 | 8.7e-11 | 2.0e-09 |
| 0.001 | 1, 1 | 1.1827579903 | 1.39891646 | 0.84548 | 6.9e-11 | 2.1e-09 |
| 0.01 | 1, 1 | 1.1898779630 | 1.41580957 | 0.84042 | 4.0e-11 | 5.2e-09 |
| 0.03 | 1, 1 | 1.2058666020 | 1.45411426 | 0.82928 | 4.2e-11 | 2.4e-09 |
| 0.1 | 1, 1 | 1.2643353456 | 1.59854387 | 0.79093 | 1.1e-10 | 4.2e-10 |
| 0.2 | 1, 1 | 1.3590805814 | 1.84710003 | 0.73579 | 5.3e-11 | 8.2e-10 |
| 0.5 | 1, 1 | 1.7097870029 | 2.92337160 | 0.58487 | 5.7e-11 | 9.5e-10 |
| 1 | 1, 1 | 2.3784104135 | 5.65683610 | 0.42045 | 1.1e-11 | 2.8e-10 |
| 2 | 1, 1 | 3.7814586005 | 14.29942915 | 0.26445 | 1.0e-10 | 4.0e-10 |
| 3 | 1, 1 | 5.2037281506 | 27.07878667 | 0.19217 | 1.1e-10 | 2.0e-09 |
| 5 | 1, 1 | 8.0613088321 | 64.98470009 | 0.12405 | 1.4e-10 | 4.9e-10 |
| 10 | 1, 1 | 15.2213845506 | 231.69054764 | 0.06570 | 1.6e-10 | 1.1e-10 |

`W² − W²_Ellis` at c = 1e-4, 1e-3, 1e-2 is 1.862863e-4, 1.864029e-3, 1.875713e-2: a quadratic fit through them has intercept 3.1e-10 and slope 1.8627 and curvature 1.298 (a three-point fit).
The rate rises with c (the AdS curvature radius is `1/√(2c)` in units of b); at c = 10 the e-folding time is 0.066 b.

Spectral check (`polar.py spectral`), `|ΔW|/W` against the shooting values at N = 240, X = 40 (S-A) and N = 360, X = 60 (S-B), the same counts at every c and no complex eigenvalue:

| c | 1e-4 | 1e-3 | 0.01 | 0.03 | 0.1 | 0.2 | 0.5 | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S-A | 6.2e-2 | 2.4e-2 | 5.7e-4 | 3.5e-5 | 5.3e-7 | 4.5e-7 | 3.9e-9 | 6.3e-11 | 2.3e-8 | 4.9e-7 | 1.3e-6 | 2.9e-4 |
| S-B | 5.9e-2 | 1.0e-2 | 6.8e-5 | 1.7e-6 | 9.0e-9 | 5.5e-10 | 4.0e-11 | 6.2e-11 | 3.4e-10 | 2.8e-9 | 6.9e-8 | **5.9e-4** |

The tolerances were declared as 2e-3 (S-A) and 2e-4 (S-B) for c ≥ 0.01; for c < 0.01 only the count is compared (the AdS end is at |ρ| > 70 and the finite domain does not resolve it). S-B at c = 10 is the one miss, reported as
not met; the count and the rate (to 6e-4) are unaffected.

**Controls** (`polar_output.txt`, `polar_derive_output.txt`; the development outcomes are in the raw README). `polar_derive.py`: `bg_V`, `lin_sign`, `ang_L`, `alg_coef`, `ellis_sign`, `pot`, `print47`, `l1_gauge`,
`susy`, `susy_end`, `code` fire on D1, D2, D4a, D3c, D4a, D5b, D5c, D6c, D6e, D6f, D7. `polar.py`: `L_sign` (23 growing modes at c = −0.3 against 0), `hz_out` (3 against 0 at c = −1.5), `rmax` (6 against 0, ℓ = 1, c = −0.3),
`g1_wrong` (3 against 0, ℓ = 1, c = −1.5), `scan_range` (0 against 1, M-AdS c = 0.5), `ads_bc` (fires the spectral check), `qnm_L` (P2), `l1_pot` (P3a, P3b), `spec_sign` (S6). Three were silent or mis-targeted in development
and were replaced or retargeted before being declared (`ang_L`, `g1_wrong`, `ads_bc`); the raw README lists them.

## 5. Where the printed text and the derivation disagree

Reported, not corrected. [ABCKY23] is arXiv:2212.12601 (the version of 27 Dec 2022; equation numbers from its PDF).

| printed | derived (this note) | what decides it |
|---|---|---|
| [ABCKY23] (47): `u′ = −(1/4r₀)[(2ω²(r²+r₀²) + (ℓ+2)(ℓ−1))K + (ℓ+2)(ℓ−1)H0 − 2iωrH1] − 2r u/(r₀²+r²)` | the same with `2ω²(r²+r₀²) − (ℓ+2)(ℓ−1)` in front of K (u = F/√2, ω = iW; the other three equations (44)–(46) are equal) | D4b: derived − printed = `K(ℓ+2)(ℓ−1)/2` exactly; D5a, b: the derived system decouples and gives the printed (52), (53), the printed would not |
| [ABCKY23] (69), ℓ = 1: `u″ = −(ω²e^{−2f} + 2/(r²+r₀²))u − 2r u′/(r²+r₀²)` | test scalar ℓ = 1: `+2/(r²+r₀²)` in place of `−2/(r²+r₀²)` (from their own (Box ψ)); the derived polar ℓ = 1 equation for the gauge-invariant scalar has neither: `χ″ + 2ρχ′/R² + ω²χ = 0` | D6d (symbolic: `ode/(R²χ″ + 2ρχ′ − W²R²χ) = −√2`); P3a, b (QNM 0.681367 − 0.617753i, the master equation at L = 2) |
| [ABCKY23] Sec. 4.2: "The polar l = 0 and l = 1 quasi-normal modes correspond to the scalar branches obtained previously [BCK18]" | the derived ℓ = 1 polar mode has the potential 1/R⁴, which is [BCK18]'s ℓ = 0 test-scalar potential (R″/R) and the Darboux partner of the printed master equation at L = 2; its QNM 0.681367 − 0.617753i is not [BCK18]'s ℓ = 1 branch (1.572708 − 0.529701i, the test scalar with 2/R² + 1/R⁴) | D6d–f, P3a, b; `l1_pot` control |

Read as wording. The ℓ = 1 polar perturbation is the backreacted one (the scalar dipole moves the throat), so its potential is not the test field's: 1/R⁴, not 2/R² + 1/R⁴. [BCK18]'s ℓ = 0 curve at M = 0 reads
ω_R = 0.684 from the plotted points of its figure (three digits; the quasinormal mode of 1/R⁴ is 0.6814 − 0.6178i), so the numbers are compatible at the reading accuracy of a figure. [BCK18] has no tables; the three-digit values read
from its EPS files are in `notes/literature.md`.

**What [BKZ12] did and did not treat.** Radial (ℓ = 0) polar and axial perturbations, the M-AdS wormholes by regularized potentials and time-domain profiles (their Figs 1, 2), the black universes by potentials and time-domain
profiles (Figs 5–8). Not: polar ℓ ≥ 1, a count of growing modes, rates. This note adds, for the same family, polar ℓ = 1, 2 (none growing at the members tested) and a count and rate of the M-AdS radial mode.

## 6. What `critical_ghost.md` §1.1 should now say

Proposed wording for the coordinator (the file is not touched here). The [BKZ12] row, result column, after the #212 update:

*Polar ℓ = 1 and 2 (#218, numerical evidence, real W in [0.01, 8] per b): no growing mode at any of the 20 black universes tested (c = −10 … −0.01, including c = −1) nor at the Ellis wormhole (there by positivity of
the master potentials, checked symbolically); complex growth rates not searched. M-AdS wormholes (c > 0, Dirichlet at the AdS end): exactly one growing radial mode at each of 12 values, c = 1e-4 … 10, rates 1.182 … 15.22 per b,
W² → the Ellis value as c → 0⁺.* And for the [ABCKY23] entry: *polar ℓ ≥ 2 master equation `V = ℓ(ℓ+1)/R² − 3/R⁴` reproduced; its ℓ = 1 statement and (47), (69) disagree with the derivation (`black_universe/polar.md`, Sec. 5).*

#211, entry on [BKZ12]: add that the M-AdS wormholes' instability is a single growing radial mode (rates above), and that polar ℓ = 1, 2 add no growing real mode for the black universes.

## 7. What this does not show

- Real W only. A growing mode with complex W (an oscillating growth) is not excluded for the black universes' polar systems, which are not self-adjoint in general (the scalar is a ghost). At the Ellis wormhole the polar problem
  is two copies of a Schrödinger operator with positive potential: none there. For the black universes the determinant has no sign change on [0.01, 8] and no interior minimum; that excludes real roots of odd and even multiplicity on the
  grid, not complex ones.
- ℓ = 1, 2 only; ℓ ≥ 3 not run. The static region ρ > ρ_h only; the Kantowski–Sachs region and formation from collapse are not touched. Members of this family and potential only, at the listed values of c.
- The M-AdS wormholes: radial modes only. Their polar modes (which need a gravitational AdS boundary condition) are not computed. W above 24 per b and below 0.01 are not searched.
- The ranges: W in [0.01, 8] (polar), [0.01, 24] (M-AdS). The count is not an argument-principle count and no rigorous bound.

## Reproduce

From this directory, with `~/Source/negative-mass-venv314/bin/python`:

```
polar_derive.py check; polar_derive.py controls; polar_derive.py emit
polar.py planted; polar.py count ELL c --res=R1|R2 --out=F.json; polar.py mads c --res=R1|R2 --out=F.json
polar.py report DIR; polar.py spectral DIR; polar.py controls
```

`polar_derive_output.txt` and `polar_output.txt` are these outputs of record, assembled from negative-mass-raw `I218-dh610`. The 108 count runs took 84 minutes of CPU on four processes (`nice -n 10`), the planted run 6 minutes, `polar_derive.py check` 133 s, the two controls runs about 25 minutes each; memory was not sampled (single-threaded numpy and numba).

## References

- [ABCKY23] B. Azad, J. L. Blázquez-Salcedo, X. Y. Chew, J. Kunz, D.-h. Yeom, "Polar modes and isospectrality of Ellis–Bronnikov wormholes", PRD 107, 084024 (2023), arXiv:2212.12601, doi:10.1103/PhysRevD.107.084024
  (INSPIRE and arXiv checked). Read in full (TeX and PDF).
- [BCK18] J. L. Blázquez-Salcedo, X. Y. Chew, J. Kunz, "Scalar and axial quasinormal modes of massive static phantom wormholes", PRD 98, 044035 (2018), arXiv:1806.03282, doi:10.1103/PhysRevD.98.044035
  (INSPIRE and arXiv checked). Read §V and its figure data; the rest not read.
- [BF06], [BKZ12]: as in `README.md` of this directory.

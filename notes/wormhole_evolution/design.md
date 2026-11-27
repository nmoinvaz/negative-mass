# A 1+1 evolution that can pass a wormhole throat: design (#39, step 1)

Step 1 of #39. No evolution code exists yet; this note, `derive.py`, `cost_model.py` and `bench/` are all there is. Status
words as in AGENTS.md: *checked symbolically* (sympy, `derive.py`), *numerical evidence*, *measured* (a timing, with its
machine). Units: a = 1 (throat radius), 8πG = 1, speed of light 1.

## Answer

Use **double-null coordinates, free evolution of (r, f, φᵢ), the second-order diamond scheme** of [SH02] and Kain [Kai26],
with D a parameter and any number of massless scalars of either sign.

- The throat is a regular point, so nothing special is done there.
- There is **no outer boundary**: the domain is the future of the corner sphere S, and its edges are null lines.
- Both ends are in the grid, and the horizons are sign tests on r₊ and r₋.
- Cost is small: one run to |p − p*| = 1e-10 takes 0.2 to 5 s of one thread, one bisection 0.2 to 4 min, for h = 0.01
  down to 0.0014 (§7).
- **Language: Python.** A numba kernel (62 ns per cell, measured) for production, and a numpy anti-diagonal version
  (192 ns per cell) as the independent reference. Rust measured 80 ns per cell, so it buys nothing here (§7).
- The hard part is not speed. **Truncation error seeds the unstable mode** (§5.3), so a "static" wormhole is static only
  for a time that grows as h shrinks. The design takes this into account: the controls tune p to p*ₕ first.
- The build is 8 to 10 h (§9). Open questions for the ruling are in §10.

## 1. Why double-null

| option | passes a throat | boundary conditions | horizon finders | cells for T = 27 (estimate) | main risk |
|---|---|---|---|---|---|
| polar-areal 1+1 (`phase5/spherical`) | no: a minimal sphere is a → ∞ | regular centre, one outer | a → ∞ | | cannot represent the data (CG §5.2) |
| Cauchy, proper distance l, gauge d = a + 2λc̄ [GGS09b] | yes | two outgoing-wave boundaries, kept outside the causal reach of the region of interest, so L ≳ T | on each slice, slicing dependent | 2L·T/(Δl·Δt) ≈ 6500/h² at Δt = Δl/4, times 3 Runge–Kutta stages | boundary pulses; slicing pathologies after collapse; outgoing-wave conditions do not preserve the constraints [GGS09b] |
| **double-null** [SH02, ST17, Kai26] | yes: r is a smooth function of (x⁺, x⁻) | none: the domain is the domain of dependence of Σ⁺ ∪ Σ⁻ | sign of r₊, r₋ | X²/h² ≈ 360/h², X = T/√2 | expanding branch may reach a finite coordinate limit (§5.5) |

The Cauchy route is a valid second code and a possible independent check of the growth rate (G1 below), not the first one.
[GGS09b] also show that the "throat as the minimum of r on each slice" is not an invariant surface, which §6 takes into
account.

[SH02], [ST17] and [Kai26] all use the dual-null form with a diamond (predictor-corrector) update and report the
static-wormhole truncation problem of §5.3. [ST17] covers n = 4, 5, 6 and Einstein–Gauss–Bonnet. Read here as data, no
code run.

## 2. Model and equations

D = n + 2. Metric ds² = −2e^{−f}dx⁺dx⁻ + r²dΩₙ². Matter: scalars φᵢ with signs εᵢ (+1 normal, −1 ghost),
G_ab = Σεᵢ(∂_aφᵢ∂_bφᵢ − ½g_ab(∂φᵢ)²), which is R_ab = Σεᵢ∂_aφᵢ∂_bφᵢ in every D (`derive.py` A2). Subscripts +, − are
∂/∂x⁺, ∂/∂x⁻.

The Ricci tensor of the full D = 4, 5, 6, 7 metric (n = 2…5 sphere angles, brute force) equals the reduced expressions in
n (`derive.py` A1). The field equations, with n symbolic (A3), are:

```
r₊₋ = −(n−1)(r₊r₋ + e^{−f}/2)/r                              (E1)
f₊₋ = −n(n−1)(r₊r₋ + e^{−f}/2)/r² + Σ εᵢ φᵢ₊ φᵢ₋           (E2)
φᵢ₊₋ = −n(r₊φᵢ₋ + r₋φᵢ₊)/(2r)                               (E3, = □φ = 0, B1)
C₊₊ ≡ r₊₊ + f₊r₊ + (r/n)Σεᵢφᵢ₊² = 0,   C₋₋ likewise          (the constraints)
```

(E1) to (E3) are evolved; the constraints are monitored. For n = 2 they are [SH02] (2.4)–(2.8).

**Identities** (*checked symbolically*: A4, B2, B3, C1, C2):
- Constraint propagation: ∂₋C₊₊ = −(n−1)(r₋/r)C₊₊ and ∂₊C₋₋ = −(n−1)(r₊/r)C₋₋, so r^{n−1}C₊₊ is constant along x⁻. A
  violation made on Σ⁺ is carried along x⁻ scaled by r^{−(n−1)}: it falls where r grows and rises only where r falls.
- Misner–Sharp mass M = r^{n−1}(1 + 2e^f r₊r₋)/2, i.e. (∇r)² = 1 − 2M/r^{n−1}. Flux laws, with the constraints
  satisfied: ∂₊M = −(rⁿ/n)e^f r₋ Σεᵢφᵢ₊² and ∂₋M = −(rⁿ/n)e^f r₊ Σεᵢφᵢ₋². M is geometric, so it does not depend on the
  normalization of φ. In D = 4, M = GM_ADM and r_AH = 2M.
- Stress components: T₊₊ = Σεφ₊², T₋₋ = Σεφ₋², T₊₋ = 0; ρ = e^f(T₊₊ + T₋₋)/2 = p_radial, and the angular pressure is
  e^f Σεφ₊φ₋. A homogeneous ghost has w = +1 with ρ < 0. So "w = −1" of #213 needs a definition of the pressure (#213
  declares it); the code outputs T₊₊, T₋₋ and the angular term.
- Expansions along the null vectors e^f∂₊ and e^f∂₋: θ₊ = n e^f r₊/r, θ₋ = n e^f r₋/r. Only their signs are used.

**Gauge.** x^± are fixed on the initial rays by the *ray gauge* f|_Σ (default f = 0). The static symmetric massless
wormhole has g_tt = −1 for every D, so with x^± = (t ± l)/√2 it is r = r(l), φ = φ(l), f ≡ 0, l = (x⁺ − x⁻)/√2 (D2, D3). Then
t = (x⁺ + x⁻)/√2 is the static time and the proper time at the throat. [SH02]'s lifetime law reads 0.60 per ln in x⁻, which is
0.846/√2: **a lifetime slope in x⁻ is 1/(√2W)**. The code reports t and the throat's proper time τ = ∫e^{−f/2}dt, never x⁻
alone. The remaining freedom x^± → F^±(x^±) changes f by ln(F₊′F₋′) and leaves the equations unchanged.

**Static solutions and exact solutions** (*checked symbolically*):

| id | solution | where used |
|---|---|---|
| D2, D3 | Ellis, any D: r′² = 1 − (a/r)^{2(n−1)}, r″ = (n−1)a^{2n−2}/r^{2n−1}, φ′ = q/rⁿ, q² = n(n−1)a^{2n−2}, ghost sign, f = 0, M = a^{2n−2}/(2r^{n−1}) → 0 at both ends | static data, controls |
| D4 | Tangherlini, any D, M ≠ 0, e^{−f} = 1 − 2M/r^{n−1}, x^± = (t ± r*)/√2 | planted check of the gravity sector with M ≠ 0, no scalar |
| D5 | **Roberts ghost solution, n = 2, f = 0:** r² = a(x⁺² + x⁻²) − x⁺x⁻, φ = √2 atan(x⁻ sinϑ/(x⁺ − x⁻cosϑ)), cosϑ = 1/(2a), a > ½. Exact and nonlinear, ghost sign. r₊ > 0 iff x⁻ < 2a x⁺, r₋ > 0 iff x⁺ < 2a x⁻, so the wedge x⁺/(2a) < x⁻ < 2a x⁺ is anti-trapped and the two lines are the trapping horizons | planted check of the full nonlinear system and of the finders, D = 4 only |
| D1 | static ODEs with a lapse ν: (e^ν rⁿν′)′ = 0 = (e^ν rⁿφ′)′, the asymmetric family of [GGS09a] | optional data, #42 |

D5 is not in the brief. It is the only exact dynamical solution with the ghost sign that I found, it crosses both r₊ = 0
and r₋ = 0 and enters the anti-trapped region, which is the structure #213 needs. A search of the literature for it was
not made. I propose it as a third planted case (§10, Q3).

## 3. The domain

```
  x⁻
   X ┌─────────────────────────────┐
     │  left end: edge x⁻ = X      │
     │                          ╱  │
     Σ⁻  (l < 0)           ╱       │ edge x⁺ = X
     │                  ╱ throat   │ (right end, l > 0)
     │               ╱  x⁺ = x⁻    │
   0 S ─────── Σ⁺ (l > 0) ─────────┘
     0                             X   x⁺
         pulse on Σ⁺ at x⁺ = x_c, moving to smaller l
```

See `domain.svg`. Data on Σ⁺ = {x⁻ = 0} and Σ⁻ = {x⁺ = 0} and at S; the evolved region is the rectangle
[0, X]². The corner S sits at the throat of the static wormhole, with r₊ = r₋ = 0 there. Σ⁺ lies in the right universe
(l = x⁺/√2) and Σ⁻ in the left one. A pulse on Σ⁺ at x⁺ = x_c is an ingoing wave and reaches the throat at x⁻ = x_c,
t = √2 x_c. The throat worldline leaves the rectangle at t = √2 X, so **the run time is set by the length of the rays**:
T = √2 X. The area of the rectangle in (t, l) equals X², since dx⁺dx⁻ = dl dt, so no cells are spent on regions the throat
cannot influence.

## 4. Initial data

All data are the characteristic data on Σ^±: f|_Σ (gauge), φᵢ|_Σ (free), r|_Σ (from the constraint), and the corner.

**Constraint solve.** On Σ⁺ the constraint C₊₊ = 0 with f₊, φᵢ₊ given is *linear* in r:

r₊₊ + f₊r₊ + (S/n) r = 0, S = Σεᵢφᵢ₊², r(0) = a, r₊(0) = 0 (the throat at the corner),

integrated with RK4 on the grid and with φᵢ evaluated analytically at the half nodes, so its error is far below the
scheme's. Σ⁻ is the same with C₋₋ and its own pulses. The corner values are f(0,0) = 0 and r(0,0) = a. The scheme needs
only grid values of r, f, φᵢ on the rays; the first derivatives are differences of them.

**Families.**
- The static wormhole: φ_static from the Ellis solution (any D), f_Σ = 0. The r that the constraint solve returns equals
  the Ellis r(l) to the integration error, because r(l) solves the same linear ODE (D3b).
- A Gaussian pulse of ghost or normal field on either ray: φᵢ(x⁺, 0) = φ_static + Aᵢ exp(−((x⁺ − x_c)/w)²), x_c ≥ 3w so
  that r₊(0) = 0 stays consistent. The pulse is in φ itself, not in ℘ = r∂φ as in [SH02] (3.1), so that it
  does not depend on D. The amplitudes are therefore not [SH02]'s c_a; the code reports each run's E₀ for comparison.
- The pulse energy E₀ is M(X, 0) − M_static(X, 0) at the far end of Σ⁺: the change of the Misner–Sharp mass along the ray.
  For a ghost pulse on the static wormhole E₀ is linear in A, for a normal pulse quadratic ([CG] §1.3).
- Optional: a two-sided symmetric pulse (the throat stays on x⁺ = x⁻ and the cost halves by reflection); Tangherlini data
  in the tortoise gauge (f_Σ ≠ 0); the asymmetric static family from the lapse ODEs.

A one-parameter family that passes through the static wormhole, such as the ghost amplitude A alone, has p* = 0 exactly
in the continuum. A "generic" family for #40 needs a second parameter (a fixed ghost pulse with E₀ < 0 and a normal pulse B swept,
or a pulse position or width swept). The code takes any parameter vector; the choice is #40's.

## 5. Scheme

### 5.1 The cell update

Cell corners S = (i, j), E = (i+1, j), W = (i, j+1), N = (i+1, j+1), Δx⁺ = Δx⁻ = h. For the vector U = (r, f, φᵢ) and
U₊₋ = F(U, U₊, U₋) from (E1)–(E3):

```
predictor:  U_c = (E+W)/2,  U₊ = (E−S)/h,  U₋ = (W−S)/h,       N̂ = E + W − S + h² F(U_c, U₊, U₋)
corrector:  U_c = (S+E+W+N̂)/4,  U₊ = (E+N̂−S−W)/(2h),  U₋ = (W+N̂−S−E)/(2h),   N = E + W − S + h² F(...)
```

One predictor and one corrector, as in [Kai26]. Order: the predictor's derivatives are first order, so N̂ has an error
O(h³); the corrector divides it by 2h, so F has an error O(h²) and N a local error O(h⁴). The (X/h)² cells accumulate
that to a global O(h²). This is the claim; §5.2 measures it. A second corrector is a parameter, to be switched on if the
measured order is below 2. F has one exponential e^{−f} and a few divisions per call; the stand-in in `bench/` has the
same arithmetic. The scheme runs by rows of x⁻, each needing only the row below, so memory is O(X/h).

### 5.2 Convergence is measured, not assumed

- Grids h, h/2, h/4 share nodes. **Self-convergence** Q = ‖U_h − U_{h/2}‖ / ‖U_{h/2} − U_{h/4}‖ → 4 for second order, over the
  region |l| ≤ L, t ≤ t_c declared per check, in max norm and L² norm.
- Against an exact solution (Ellis static, Tangherlini, Roberts): the error itself, ratio 4 per halving.
- **Constraints:** C₊₊ and C₋₋ at cell centres (second-order differences of the stored values) must fall by 4 per halving.
- **Mass-flux identity:** M(X, x⁻ₐ) − M(X, x⁻_b) against the integral of the flux along the edge, same order.
- Near threshold: p*ₕ at three h, and the slope of lifetime against ln|p − p*ₕ| at three h. The order of p*ₕ is whatever
  it measures.

### 5.3 Truncation seeds the unstable mode

The discrete scheme is a map with a fixed point (the discrete static wormhole) that differs from the sampled exact static
wormhole by O(h²). That difference has a component along the unstable mode. So an unperturbed run leaves the wormhole after
t ≈ (1/W) ln(1/(c h²)): [SH02], [ST17] and [Kai26] all observe it (SH02 Fig. 2 needs up to 9601 points on x⁺ ∈ [0, 20]
to keep it static for the time shown). If c ≈ 1 (not measured yet), h = 0.01 at D = 4 gives about 7.8 a; c is measured in step 2.

Consequences, built into the design:
1. **"Static to the truncation error" is a windowed statement:** the deviation of an unperturbed run converges at
   second order at each fixed t, for t up to a declared t_c (≈ 4–5 a, where e^{t/T} ≤ 10²). Beyond t_c the deviation grows
   as e^{t/T}, which is physics, and it is itself a measurement of c.
2. **Bisection absorbs it.** In a family p the data are the sum of the family's perturbation and this fixed O(h²) seed.
   The threshold of the *discrete* system is p*ₕ = p* + O(h²). The lifetime law keeps the slope 1/W up to O(h²) in the
   discrete growth rate. #40 and #150 report p*ₕ at three h and the slope at three h.
3. **The growth-rate control (G1) tunes first.** A single run at small amplitude cannot measure W: the seed puts a floor
   of ~1e-4 under the unstable amplitude. G1 bisects to p*ₕ with |p − p*ₕ| ≤ 1e-9, then fits the departure
   ln|r_th − a| against τ over a window in which the amplitude runs from 1e-6 to 1e-3 (window to be declared after a
   development run). This is a clean linear growth over 7 e-folds.
4. Roundoff seeds ~1e-16, six decades below |p − p*| = 1e-10, so the bisection to 1e-10 is not limited by it.

### 5.4 Horizons, excision, stop conditions

- Trapped: r₊ < 0 and r₋ < 0. Anti-trapped: r₊ > 0 and r₋ > 0. Untrapped right of the throat: r₊ > 0 > r₋, left of it
  the reverse. The throat sphere of the static wormhole has r₊ = r₋ = 0.
- Contours r₊ = 0 and r₋ = 0 by linear interpolation on cell centres, each with the tangent dx⁺/dx⁻: > 0 timelike,
  < 0 spacelike, 0 or ∞ null (a trapping horizon's type, Hayward's classification as used by [SH02]).
- **Excision** as in [SH02]: a cell with r < r_cut or e^{−f} > F_cut (defaults 0.05 a and 25, to be declared before the
  first run of record) is dead, and so is every cell with a dead S, E or W. A dead cell inside the apparent horizon cannot
  influence anything outside the event horizon, so there the excision is exact. A cell that e^{−f} > F_cut kills outside
  the event horizon, where the grid has stretched, truncates the domain: nothing in its future is reported. The first
  dead cell of each run and whether it was inside a trapped region are recorded.
- A run **stops** when its fate is classified (declared per issue): trapped region present and its outer boundary
  r₊ = 0 null or spacelike for a declared length; or anti-trapped region present with r_th ≥ 2a; or the rectangle is
  exhausted. A run that reaches the rectangle's end without either is reported as "static through X".

### 5.5 Where this form can fail (stated before it is tried)

- **Grid stretching** at a horizon: after the horizon forms, e^{−f} grows and resolution in x⁻ falls. The mass at
  *formation* is clean. A mass after the horizon has settled needs a stated tolerance, and is limited by this
  (and by the excision).
- **The expanding branch** may reach a finite coordinate limit: in de Sitter-like expansion future infinity lies at finite
  null coordinates. The code records max e^{−f} and the e-folds of r_th reached per run, so the limit is measured. A
  gauge on the rays that compresses far coordinates (F^±) is the first remedy. I cannot say now how many e-folds fit.
  [GGS09b] followed an exponential expansion in their Cauchy gauge, and [SH02] report H ≈ 1.1/a in this form.
- **Regular centres** (the polar-areal data of #45 and #46) are not in this domain. A regular centre needs a boundary line
  r = 0 with its own treatment. It is an addition to the same cell update, not part of the minimal build.

## 6. What #40, #150 and #213 get

| need | the code outputs | how it is checked |
|---|---|---|
| lifetime, horizon time | per run: t and x⁻ of the first trapped cell, of r_th = a(1 ± 0.1…2), of the departure; the contours | Roberts (D5) places r₊ = 0, r₋ = 0 exactly: finder error O(h²) and a mutated finder fails |
| trapped, anti-trapped regions, horizon type | masks per row, contour list with tangent, H′ | D5 anti-trapped wedge with known edges |
| mass at each end | M(X, x⁻) on edge x⁺ = X (right) and M(x⁺, X) on x⁻ = X (left), minus M_static; [Bondi-type: M at the outermost point of each outgoing ray, minus the static tail a^{2n−2}/(2r^{n−1}). The tail is 0.036 a at r = 14 in D = 4, 12 % of [SH02]'s threshold mass 0.30 a, so it is subtracted, not neglected] | flux identity (§5.2); static Ellis M tends to 0 |
| throat radius and expansion rate | on the locus r₊ = r₋ (the minimum of r on the slices t = const of the ray gauge): r_th(t), τ, H = e^{f/2}(r₊ + r₋)/(√2 r) with u = e^{f/2}(∂₊ + ∂₋)/√2; also the value on the diagonal x⁺ = x⁻ [SH02]. The two differ for one-sided data; [GGS09b] show the minimum is slicing dependent, so both are written | static: r_th = a, H = 0; D5 has a closed form |
| ρ, p, w | T₊₊, T₋₋, angular term on the throat locus | formulas of §2 (checked symbolically); definition of w is #213's |
| black-hole mass | M at the trapped region's outer edge at formation: M = r^{n−1}/2 on r₊ = 0 | Tangherlini planted: horizon at r = (2M)^{1/(n−1)} |
| a run long enough | X is a parameter; per-run max e^{−f}, e-folds, and the first excised cell | cost table, §7 |
| a normal field beside the ghost | εᵢ list, any number of fields; the normal field decouples from the static background | #213's control: a positive-energy pulse collapses |

The code streams the diagnostics row by row and stores the two edges, the throat series, the contours and a decimated grid,
so memory does not grow with the grid.

## 7. Cost

The kernel stand-in (`bench/`) has the arithmetic of §5.1 for four fields, and is run on exact Ellis data at h = 0.001
where the static solution has not begun to leave. It is not the evolution code and checks nothing. dh610-fedora,
i7-12700, load average 9 of 20 (shared with choptuik-large-d), nice 10, **one thread**, best of three at N = 2000:

| implementation | ns per cell | Mcell/s |
|---|---|---|
| numba (row loop) | 61.7 | 16.2 |
| Rust (row loop, `-C target-cpu=native`) | 79.6 | 12.6 |
| numpy (anti-diagonals) | 192.3 | 5.2 |

numpy is 3× slower than numba because each of its anti-diagonal sweeps is a few thousand elements long, so the interpreter
overhead per array operation is amortized only partly. The row loops of numba and Rust are limited by the latency of the
W → N dependence in the row (an anti-diagonal Rust loop could be faster; not tried).

`cost_model.py` (output `cost_model_output.txt`): run length T = √2·3 + ln(1/|p − p*|)/W + 3, with |p − p*| = 1e-10, and
45 runs per bisection (about 33 halvings plus the lifetime points; no early exit counted, which would save about a third).
For h at D = 4 and scaled by 1/√(D − 3) at larger D:

| D | 1/W | T | X | h | cells | numba per run | numpy per run | numba per bisection | numpy per bisection |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 0.846 | 26.7 | 18.9 | 0.01 | 3.6e6 | 0.22 s | 0.69 s | 0.2 min | 0.5 min |
| 4 | 0.846 | 26.7 | 18.9 | 0.0025 | 5.7e7 | 3.5 s | 11 s | 2.6 min | 8.2 min |
| 5 | 0.579 | 20.6 | 14.5 | 0.0071 | 4.2e6 | 0.26 s | 0.81 s | 0.2 min | 0.6 min |
| 5 | 0.579 | 20.6 | 14.5 | 0.0018 | 6.8e7 | 4.2 s | 13 s | 3.1 min | 9.8 min |
| 6 | 0.462 | 17.9 | 12.6 | 0.0058 | 4.8e6 | 0.30 s | 0.92 s | 0.2 min | 0.7 min |
| 6 | 0.462 | 17.9 | 12.6 | 0.0014 | 7.7e7 | 4.7 s | 15 s | 3.6 min | 11 min |

One bisection per (D, h) over all 12 (D, h) pairs of the full table costs 0.21 h of one thread in numba and 0.64 h in
numpy, per family. At 4 threads the bisections are independent, so the wall time is a quarter. **The bisections of #40 and #150 (two
families, three D, three resolutions) fit in under an hour of 4 threads in either language.** The choice of language is
therefore about reliability and effort.

**Recommendation: Python.** It keeps the sympy-derived formulas, the checks, the `--mutate=` hooks and the analysis in one
language, and numba==0.68.1 is already pinned in `requirements.txt` for `phase5/spherical/evolve.py`. The numba kernel is
production; the numpy anti-diagonal version is a second implementation of the same scheme used in the planted checks, so
an indexing error in one shows as a difference between the two (it does not test the scheme itself, which the exact
solutions do). **Rust** would be the choice if a runner cannot install numba (rpi5, Windows are runners in this project);
its kernel is 80 ns per cell. I have not written either.

## 8. Step 2: checks and their controls (proposed; the README of the raw folder gets the final wording after the ruling)

Every check can fail, and each has a `--mutate=` whose targets are computed in a development run before they are declared
(AGENTS.md). Targets for `derive.py`'s mutations were computed that way: dev_05 to dev_07 in the raw folder, and all nine fire exactly
their declared checks in run02.

| id | check | planted truth / control | proposed pass | mutation (declared after dev run) |
|---|---|---|---|---|
| P1 | static Ellis stays static to the truncation error | exact, D = 4, 5 | deviation at t ≤ t_c converges at order 2 ± 0.25 at h, h/2, h/4 | sign of the sphere term; ghost sign |
| P2 | Tangherlini stays Tangherlini, M ≠ 0, D = 4, 5, 6, exterior only | exact, vacuum, no unstable mode, so the order is measured over long times with no seeding (§5.3) | same order | sphere coefficient |
| P3 | Roberts ghost solution, D = 4 (needs ruling Q3) | exact, nonlinear, ghost | error order 2; finders return the lines x⁻ = 2a x⁺, x⁺ = 2a x⁻ to O(h²) | scalar coefficient; finder sign |
| C1 | C₊₊, C₋₋ converge to 0 at order 2 | by (B2) | ratio 4 ± 1 | constraint sign |
| C2 | mass flux identity at order 2 | by (C1, C2) | ratio 4 ± 1 | mass exponent |
| G1 | perturbed wormhole grows at 1/W: tune p to p*ₕ, fit the departure | W = 1.181970 (D = 4), 1.727703 (D = 5) from #195's W²a² | within 1 % | sphere coefficient; time unit √2 |
| F1 | fates at D = 4: one pulse of each sign | [SH02]: E₀ < 0 expands, E₀ > 0 collapses, a normal pulse collapses. Which sign of A has E₀ < 0 is read from the run's E₀, not assumed | qualitative, reported as seen | none (not tuned) |

The controls G1 are the two the issue asks for. Their window and fit rule are declared after a development run, together
with the thresholds of §5.4, and then not moved.

## 9. Estimate

| part | hours |
|---|---|
| kernel (numba and numpy), initial data (static ODE, constraint solve, pulses, Tangherlini, Roberts) | 3 |
| diagnostics: finders, throat, mass, flux, constraints, streaming output | 2 |
| checks P1–P3, C1–C2, G1, F1 with `--mutate=`, development runs, then runs of record | 3.5 |
| note, raw folder, posts | 1.5 |
| **total** | **8 to 10** |

Step 1 took about 1 h of wall time against the 2 to 4 h allowed. The kernel stand-in is not reused.

## 10. Questions for the ruling

1. **Formulation:** double-null, free evolution, as in §1 to §5. A Cauchy proper-distance code is the alternative.
2. **Language:** Python (numba production kernel plus numpy reference) as in §7, or Rust.
3. **Extra planted cases** beyond the brief: Tangherlini (P2) and the Roberts ghost solution (P3). Both are cheap. P3 is the
   only check of the nonlinear ghost dynamics and of the finders against an exact answer. Without P3 the only nonlinear
   check is the qualitative F1.
4. **Throat definition for H:** write both the minimum of r on the ray-gauge slices and the diagonal x⁺ = x⁻, and let #213
   choose (§6). Or fix one now.
5. **Excision and stop thresholds** (r_cut, F_cut, the "classified" rules of §5.4): declare in the README before the first
   run of record, as proposed, or fix them here.
6. **G1 protocol:** tune p to p*ₕ then fit, because a single small-amplitude run is at the mercy of the O(h²) seed (§5.3).
   If that does not count as "a perturbed wormhole grows at the linear rate", say so before the build.
7. **Regular centre** for #45 and #46: not in scope here (§5.5). Leave it, or add the centre boundary to this build.

## References (read here as data, none run)

- [SH02] H. Shinkai, S. A. Hayward, PRD 66, 044005 (2002), arXiv:gr-qc/0205041. Dual-null form (2.3)–(2.17), grid and
  scheme (II B), static data (2.34)–(2.35), pulses (3.1), mass (2.36)–(2.39), excision at 25.
- [ST17] H. Shinkai, T. Torii, PRD 96, 044009 (2017), arXiv:1706.02070. n-dimensional dual-null with Gauss–Bonnet; the same
  grid, the same excision, and the same remark that truncation error destroys the static wormhole.
- [GGS09a] J. A. González, F. S. Guzmán, O. Sarbach, CQG 26, 015010 (2009), arXiv:0806.0608. The linear mode (see
  `ellis_efold.md`).
- [GGS09b] J. A. González, F. S. Guzmán, O. Sarbach, CQG 26, 015011 (2009), arXiv:0806.1370. Cauchy form (1)–(9),
  gauge d = a + 2λc̄, radiative boundaries pushed outside the light cone, second-order convergence, the throat as a slicing
  dependent minimum.
- [Kai26] B. Kain, arXiv:2602.18231. A pedagogical dual-null code: predictor-corrector diamond update (its §2.3), grid
  spacing Δu = Δv = 1/N, static configuration lasting longer as the spacing falls.
- [CG] `notes/critical_ghost.md` §1.2, §1.3, §5.2–5.3. #195: `notes/wormhole_largeD/w2_output.txt` (W²a²).

Files: `derive.py` (25 checks; 9 mutations; log `derive_output.txt`; raw `I39-dh610/derive/run02_*`), `cost_model.py`,
`bench/` (stand-in kernels, raw `I39-dh610/bench/`), `domain.svg`.

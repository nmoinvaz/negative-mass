# Phase 3 results — signed point particles at first post-Newtonian order

Model: S = Σ_A s_A(−m_A ∫dτ_A) + S_EH, m_A > 0, s_A = ±1. This is the point-particle version of Shirokov's
S = ∫√−g [R/16π + L(Φ₊) − L(Φ₋)]. Write **μ_A ≡ s_A m_A**. G = c = 1. All checks are in this directory, and each
script `X.py` saves its output to `X_output.txt`.

| script | what it checks |
|---|---|
| `point_particle.py` | worldline EL equations in an arbitrary metric + EM potential; p, E; Lorentz-force response |
| `eih.py` | EIH Lagrangian → EL → order-reduced 1PN acceleration = Blanchet eq. (203); E, P, J, G conserved (signed masses) |
| `test_body.py` | Blanchet's acceleration in the test-body limit = geodesic of harmonic-gauge Schwarzschild to 1PN |
| `circular_check.py` | ω²(r) and the gauge-invariant E(x) at 1PN against the literature, including a mixed-sign mass pair |
| `collinear.py` | collinear pair / dipole formulas; the condition for a rigid pair; Shirokov's numbers |
| `integrate_1pn.py` | numerical 1PN trajectories (exact EL and order-reduced); conservation; separation drift |
| `breakdown.py`, `radiation.py`, `overlap.py` | where 1PN breaks down; Newtonian multipole radiation; Newtonian finite-size overlap of the Phase 2 stars |

## 1. Which masses are negative

The action depends on the particles only through μ_A, and depends on it linearly.

- **(a) Passive mass / the geodesic coefficient.** The EL equations of −μ∫dτ + q∫A in an arbitrary metric reduce
  to μ u^ν∇_ν u^μ = q F^μ_ν u^ν. The symbolic residual is 0 (`point_particle.py`). With q = 0 the equation is
  μ × (geodesic equation) = 0, so **the sign of μ has no effect**. "Passive mass" is the coefficient multiplying
  ∇Φ in m_i a = −m_p ∇Φ, and it equals μ, the same as the inertial coefficient. Only the ratio is observable, and
  that ratio is +1.
- **(b) Active mass.** T^{μν} = Σ μ_A ∫u^μu^ν δ⁴/√−g dτ_A, so the metric is sourced by μ_A. The ADM/Komar mass of
  a lone particle is μ_A < 0 for s_A = −1.
- **(c) Inertial mass in the action.** L = −μ√(1−v²) gives p = μγv and E = μγ. For μ < 0 both are negative, and
  so is the kinetic term μv²/2. The conserved momentum is the signed one.
- **(d) Response to a non-gravitational force.** This depends on where the charge coupling sits.
  - **Standard-sign coupling** (q∫A outside the flipped sector): μ d(γv)/dt = q(E + v×B). A μ < 0 body accelerates
    *against* the applied force. This is Bondi's negative inertial mass.
  - **Coupling inside the flipped sector** (s_A(−m∫dτ + q∫A)), which is what Shirokov's sign-flipped Lagrangian
    gives if Φ₋ is charged: the EOM is m d(γv)/dt = q(E + v×B), so the body responds normally. However, it
    *sources* the Maxwell field with charge s_A q = −q. That is exactly a Bondi body (μ = −m) carrying charge −q,
    with the same q/μ ratio.

  Either way, **every gravitational and electromagnetic observable of a phantom point particle is that of a Bondi
  all-negative mass μ = −m.** The only free choice is which charge label you attach.

**Comparison with Bondi's cases.**

- *All three masses negative* (Bondi): this is exactly the model, with μ_A = s_A m_A.
- *Positive inertial, negative active* (Bonnor's footnote case): in this model it is **not a different theory**,
  only a different bookkeeping. Shirokov's "positive inertial and passive mass" refers to the matter field
  equation of Φ₋ being unchanged. For the worldline that equation is μ × (geodesic) = 0, which the sign of μ does
  not affect. If you call the inertial mass +m, the unsigned momentum m(v₁ + v₂) of the dipole grows (the third
  law fails). If you call it −m, the momentum m v₁ − m v₂ = 0 is conserved. The motion is the same in both
  descriptions, and only the second one is a Noether charge.
- **Is the phantom distinguishable from Bondi's all-negative mass using gravity alone? No.** The expectation is
  confirmed, not assumed: the reduced (Fokker) action is the m_A → μ_A continuation (§2), so every gravity-only
  observable is identical. A difference can only appear through non-gravitational couplings, and even there it is
  a relabelling of charges (case d).

## 2. 1PN (EIH) Lagrangian and equations of motion

**Why EIH with m_A → μ_A is correct.**

1. The field equation is G_μν = 8π Σ_A μ_A T_A,μν[z_A]. The near-zone iteration (h ~ μ/r, then (μ/r)², …) is a
   power series in the μ_A with no |μ_A| anywhere.
2. The matter equations are ∇·T_A = 0, which is the same as the geodesic equation independent of μ (§1a). The
   per-body equations are therefore those of standard EIH with m → μ.
3. Self-field divergences (terms like μ₁²μ₂/r²) are removed by Hadamard / dimensional regularization. This
   operation is polynomial in μ, so the regularized coefficients continue to μ < 0.
4. For extended phantom bodies, "effacement" of internal structure at 1PN only uses ∇·T_eff = 0 and the virial
   relation for a stationary T_eff. Both hold for T_eff = −T(Φ₋), because Φ₋ obeys its own unchanged field
   equation.

**Lagrangian** (Landau–Lifshitz eq. 106.17; Blanchet, Living Rev. Relativ. 17, 2 (2014), 1PN harmonic gauge):

    L = Σ μ_A v_A²/2 + Σ μ_A v_A⁴/8 + μ₁μ₂/r
        + μ₁μ₂/(2r) [3(v₁² + v₂²) − 7 v₁·v₂ − (n·v₁)(n·v₂)] − μ₁μ₂(μ₁ + μ₂)/(2r²)

**Acceleration** (Blanchet eq. 203 at 1PN, m → μ), where n = (x₁ − x₂)/r:

    a₁ = −μ₂ n/r² + (1/c²){ [5μ₁μ₂/r³ + 4μ₂²/r³ + (μ₂/r²)(3/2 (n·v₂)² − v₁² + 4 v₁·v₂ − 2v₂²)] n
                           + (μ₂/r²)[4 n·v₁ − 3 n·v₂](v₁ − v₂) }

**Verification.**

- **(i) EL equations vs. eq. 203.** I took the EL equations of L and order-reduced them (accelerations inside the
  O(c⁻²) terms replaced by Newtonian ones). The result agrees with eq. 203, with residual ≲ 10⁻⁴⁶ at random
  points with μ₁ > 0 > μ₂, for both the O(1) and O(c⁻²) coefficients (`eih.py`).
- **(ii) Test-body limit.** Eq. 203 at μ₁ → 0 equals the geodesic of harmonic-gauge Schwarzschild expanded to
  1PN, which I derived from the Christoffel symbols. The residual is exactly 0 (`test_body.py`).
- **(iii) The μ₁μ₂ terms.** These are invisible in (ii), so I checked them separately. Circular orbits give
  ω² = (M/r³)[1 + (ν − 3)M/r], and E(x) = −(Mνx/2)[1 − (3/4 + ν/12)x]. Both are exact, and both hold for a
  mixed-sign pair (m₁ = 2, m₂ = −1/3) (`circular_check.py`).

This verifies the EIH starting point against the literature form.

## 3. Bondi dipole: μ₁ = +m (front, at x = d), μ₂ = −m (rear, at x = 0)

**Released from rest** (exact algebra, `collinear.py`):

    a₁ = (m/d²)(1 − m/d),   a₂ = (m/d²)(1 + m/d),   a_mid = m/d²  (no 1PN correction),
    d̈ = a₁ − a₂ = −2m²/d³.

**With a common velocity V** along the axis:

    a₁,₂ = (m/d²)(1 − 5V²/2 ∓ m/d),   a_mid = (m/d²)(1 − 5V²/2).

**With a relative velocity u:** d̈ = −(m/d²)(Vu + 2m/d).

**With a common transverse velocity V:** a₁,₂ = (m/d²)(1 − V² ∓ m/d).

The −5V²/2 factor is pure kinematics. Converting to proper acceleration (×γ³) and rest-frame separation
(d₀² = d²γ²) gives (1 − 5V²/2)(1 + 3V²/2)(1 + V²) = 1 + O(V⁴). The numerical values at V = 0.05, 0.1, 0.2 are
0.99999, 0.99981 and 0.9967 (`integrate_1pn.py`). The 1PN dipole is therefore Lorentz-covariant: there is
**no 1PN correction to the proper acceleration m/d²**. The only 1PN corrections are ∓m/d on the two bodies
(±0.14 % at d = 10, m = 0.0143).

**Is the pair rigid? Not for matched masses.**

- In harmonic coordinates, d̈ = −2m²/d³ at rest. Born rigidity in a frame accelerating at a_mid requires
  d̈ = −a_mid² d = −m²/d³.
- The proper-separation acceleration is therefore d̈ + a_mid² d = −m²/d³ < 0: **the matched dipole slowly
  contracts.**
- General result:

      (proper separation)¨|_rest = −E_tot/d²,   E_tot = μ₁ + μ₂ − μ₁μ₂/d

  (exact to O(c⁻²) when μ₁ + μ₂ = O(c⁻²); `collinear.py`). E_tot is the 1PN total mass-energy, i.e. the rest
  masses plus the (positive, for opposite signs) interaction energy.
- **Rigid ⇔ E_tot = 0 ⇔ |μ₂| = μ₁(1 + μ₁/d).** The negative mass must be heavier by μ₁²/d, i.e. 0.143 % at
  d = 10.

Why E_tot = 0 is required: a boost-invariant (uniformly accelerated) system must have zero total four-momentum.
At 1PN this shows up as the centre-of-energy integral K = G − Pt being constant. A rigid runaway with E_tot ≠ 0
would make K grow like E_tot·X_mid(t).

**Numerical check** (`integrate_1pn.py`, d₀ = 10, m = 0.0143, values of d·γ = proper separation):

| case | t = 200 | 400 | 1000 | 2000 |
|---|---|---|---|---|
| matched | 9.9959 | 9.9836 | 9.8979 | 9.592 |
| \|μ₂\| = μ₁(1 + μ₁/d) (rigid) | 9.99999 | 9.99998 | 9.99985 | 9.9991 |
| \|μ₂\| = 0.9961 μ₁ | 9.985 | 9.939 | 9.616 | 8.42 |
| \|μ₂\| = 1.0039 μ₁ | 10.007 | 10.028 | 10.174 | 10.68 |

- The matched case follows δ = −m²t²/(2d³) = −0.102 at t = 1000.
- In the rigid case, the residual drift at t = 2000 (v = 0.28) is O(c⁻⁴).
- Hyperbola centres for the two bodies stay together to 1 part in 10⁶ (rigid case) but separate in the other
  cases.

**Mismatched masses.** The drift is set by the sign of E_tot:

- If the negative mass is too light (E_tot > 0), the pair closes up.
- If it is too heavy (E_tot < 0), the pair separates.
- At rest the drift is δ̈ = −E_tot/d², growing to |Δμ|/d² at leading order.

**Comparison with Shirokov.**

- **The a·d²/M̄ ≈ 1.08 at d = 10 is not explained by 1PN.**
  - The midpoint has no 1PN correction, and the individual bodies are off by only ±m/d = ±0.14 %.
  - With his quoted masses (M₊ = 0.014350, M₋ = −0.014295) I get a·d²/M̄ = 0.9966 (front) and 1.0033 (rear).
  - Newtonian finite-size overlap of the actual Phase 2 star profiles *reduces* the force, by 0.085 % at d = 10
    and 0.9 % at d = 8 (`overlap.py`; 2.5 % of each star's mass lies beyond d/2). So it has the wrong sign and is
    100× too small.
  - The Phase 2 mass-normalisation gap (our ADM masses are 2 % higher) would account for ≈ +2 % if his M̄
    underestimates the true mass. That is still not 8 %.
  - The 1PN effective exponent is −2 ± m/d = −2 ± 0.0014 per body, not −2.03.
  - What is left is a ~4 % systematic in the measured d (tracker or coordinate definition), or strong-field and
    structure effects not captured here. I could not read the paper's methods: the arXiv PDF has no extractable
    text.
- **Late-time separation growth.**
  - With his quoted masses, E_tot = +7.6 × 10⁻⁵ > 0, so **1PN predicts slow contraction**:
    δ(200) = −0.015, δ(400) = −0.060. The Phase 2 masses give δ(400) = −0.054.
  - Growth would need E_tot < 0. The quoted N = 192 value d(400) ≈ 10.7 (+0.7) would need E_tot ≈ −8.8 × 10⁻⁴,
    i.e. |M₋| larger than M₊ by ~6 %, far outside the stated 0.39 % match.
  - **1PN therefore supports his attribution of the growth to grid error.** Caveat: coordinate separation is
    gauge-dependent, and his gauge is not harmonic. At 1PN, though, a gauge change only shifts d by O(m) ≈ 0.01,
    not 0.7.

## 4. Conserved quantities and Bergmann

The signed system has the standard EIH integrals with m → μ (p_A ≡ ∂L/∂v_A):

    E = Σ μ_A v_A²/2 − μ₁μ₂/r
        + (1/c²)[ 3/8 Σ μ_A v_A⁴ + μ₁μ₂/(2r)(3v₁² + 3v₂² − 7v₁·v₂ − (n·v₁)(n·v₂)) + μ₁μ₂(μ₁ + μ₂)/(2r²) ]

    P = Σ μ_A v_A [1 + (1/c²)(v_A²/2 − μ_B/(2r))] − (1/c²) μ₁μ₂/(2r) [(n·v₁) + (n·v₂)] n

    J = Σ x_A × p_A,   G = Σ μ_A x_A [1 + (1/c²)(v_A²/2 − μ_B/(2r))],   with Ġ = P.

**Symbolic check** (`eih.py`): on-shell dE/dt, dP/dt and dJ/dt, and d(G − Pt)/dt, vanish through O(c⁻²), with
residuals ≲ 10⁻⁴⁶ at random signed-mass points. A control without the 1PN terms in G fails (residual 712).

**Numerical check** (dipole, t ≤ 2000, v up to 0.28):

| integrator | E | P | J (off-axis run) | K = G − Pt |
|---|---|---|---|---|
| exact EL of L | drift ~ 10⁻¹⁸ (machine precision) | ~ 10⁻¹⁸ | ~ 10⁻¹⁵ | O(c⁻⁴) |
| order-reduced eq. 203 | 1.3 % of m²/d | 10⁻⁴ m | O(c⁻⁴) | O(c⁻⁴) |

For the exact EL integrator, E, P and J are exact Noether charges of L. K drifts by ≈ 0.006–0.007 m·d at v = 0.3
because L is Poincaré-invariant only to O(c⁻²); this is an O(c⁻⁴) effect. For the order-reduced integrator, every
quantity is conserved up to O(c⁻⁴), consistent with truncation, and dP/m scales ∝ m/d as expected.

**The physics.**

- The dipole has E = m²/d > 0 (interaction energy only, since the rest masses cancel) and P = 0 exactly.
- The pair speeds up with no expenditure of energy or (signed) momentum: E_kin = (μ₁ + μ₂)v²/2 = 0.
- Unsigned momentum m(v₁ + v₂) grows linearly in t. That is the "third-law violation" seen by anyone who counts
  inertial mass as positive.

**Bergmann.** His claim is that active and passive mass cannot differ in GR. **It holds in this model.** Each body
has one parameter μ_A, which is simultaneously its inertial (§1c), passive (§1a) and active (§1b) mass. The model
is Bondi's "all three negative" body, not a body with active ≠ passive.

- **What Bergmann's argument does not exclude:** the runaway itself, and the *sign* of all three masses being
  negative.
- **Where an apparent violation comes from:** Shirokov's "positive inertial/passive, negative active" description
  is a statement about the matter field equation. It does not carry over to an effective point mass. The
  apparent active ≠ passive arises only from labelling the inertial mass by the sign of the Φ₋ kinetic term in
  its own equation, which is invisible to gravity.
- **Third law:** it holds in signed form (dP/dt = 0 with P = Σμv) through 1PN, and Bonnor's "violation" is the
  unsigned bookkeeping.
- **Caveat:** if a non-gravitational force with *standard-sign* coupling acts on the phantom (case d, first
  sub-case), the body responds as a negative inertial mass. That is still Bondi, not Bergmann-violating.

## 5. Validity and breakdown

- **Compactness.** 1PN needs m/d ≪ 1 (here 1.4 × 10⁻³) and v² ≪ 1.
- **Speed.** The runaway speed grows without bound, approaching c in proper time. In coordinate time it follows
  hyperbolic motion v = αt/√(1 + α²t²) with α = (m/d²)(1 + m/2d) for the rigid pair.

1PN integration of the rigid pair compared with exact hyperbolic motion (`breakdown.py`):

| t | 1000 | 2000 | 3000 | 4000 | 5000 | 6000 | 7000 |
|---|---|---|---|---|---|---|---|
| v (hyperbola) | 0.142 | 0.275 | 0.394 | 0.497 | 0.582 | 0.651 | 0.708 |
| rel. error of v_1PN | 1.5e-4 | 2.1e-3 | 8.6e-3 | 2.1 % | 4.1 % | 6.7 % | 9.9 % |

- The error scales as v⁴ (2PN), going from ×14 between v = 0.14 and v = 0.28.
- Proper separation, which is constant in exact rigid motion, degrades from 9.9999 to 9.95 at t = 5000 and to
  9.09 at t = 8000. The 1PN velocity eventually overshoots past v = 0.85.
- **Estimate: 1PN is good to ~1 % for v ≲ 0.4 and fails (≳ 10 %) for v ≳ 0.7, i.e. t ≳ 0.4–0.7 d²/m (≈ 3000–5000
  for Shirokov's parameters).** His runs (t ≲ 400–1000, v ≲ 0.06–0.14) are well inside the 1PN regime.
- This test treats the boost-symmetric exact solution as hyperbolic motion. That assumption is motivated by
  Bonnor–Swaminarayan but not proven here for this system.

**Radiation** (`radiation.py`, Newtonian multipoles about a fixed origin):

- The mass dipole Σμx = m d is constant, and the quadrupole has a vanishing third derivative.
- The first radiating moment is the octupole: L₃ = 24 m⁶/(35 d⁶) = 5.9 × 10⁻¹⁸ at d = 10. That is ~10⁻¹⁴ of E_tot
  per 1000 time units, consistent with "no detectable radiation".
- Speculative: losing energy would push E_tot < 0 and drive slow separation *growth* (δ̈ = −E_tot/d²), but far too
  little to matter.
- The fixed-origin slow-motion expansion itself fails once the displacement is comparable to a wavelength
  (v → 1), so this is an order-of-magnitude statement only.

## Caveats / not verified

- I did not read Bondi 1957, Bergmann's discussion or Shirokov's methods directly. The Shirokov PDF returned no
  text. The comparison in §3 rests on the numbers quoted in the task and in `notes/literature.md`.
- The EIH formula is checked against Blanchet eq. 203 as I transcribed it, and *independently* by (i) the
  geodesic limit and (ii) the gauge-invariant E(x) and ω²(r). The transcription is therefore confirmed even though
  the Living Review text was not fetched.
- Shirokov's acceleration and separation are measured in his gauge. 1PN predictions are in harmonic gauge.
  Gauge-invariant statements here are a_proper, proper separation, E_tot and P.
- 1PN, including the "rigid ⇔ E_tot = 0" result, is perturbative. The statement that a boost-symmetric system
  needs zero four-momentum is general, but I have not proven that the exact GR dipole with E_tot ≠ 0 drifts as
  −E_tot/d² beyond 1PN.

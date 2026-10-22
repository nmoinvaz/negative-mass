# Phase 4 addendum: gravitational radiation from the runaway ± pair

**Question.** Can gravitational radiation drive the Bondi mass of a self-accelerating ± pair to −∞? Phase 4
proved only the per-slice bound M_B(u) ≥ −W₋(u).

**Answer in brief.**

1. **An eternally uniformly accelerated pair radiates without limit.** In linearized gravity the radiated
   4-momentum is dP_rad^μ/dτ = L₀ u^μ. Here L₀ = (24/35) m⁶/d⁶ (= 5.9×10⁻¹⁸ for m = 0.0143, d = 10), with u^μ
   the pair's 4-velocity. This is the gravitational analogue of Larmor's result for a hyperbolic charge. In a
   fixed Bondi frame the Bondi energy falls like E_B(t) = E₀ − L₀ t. That is unbounded, but very slow:
   2.4×10¹⁵ time units to lose one m.
2. **The exact solution agrees.** The Bonnor–Swaminarayan (BS) strut-free solution is the exact boost-rotation-
   symmetric version of this pair. Its news has exactly the linearized angular and time dependence (up to an
   overall factor −2 that I have not resolved). On the regular half of 𝒥⁺ (u > 0), its Bondi mass is finite,
   positive and decreasing, with m(u) → 0 as u → ∞ and m ~ (3/8)D²A²/u → +∞ as u → 0⁺. The infinite radiated
   energy piles up at the corner (u = 0, θ = 0) of 𝒥⁺, where the particle worldlines end.
3. **No contradiction with Phase 4.** The floor −W₋ ≈ −e^η m recedes at the same exponential rate in rapidity
   as the radiated energy L₀t ≈ L₀e^η/(2A). Their ratio tends to (12/35)(m/d)⁴ = 1.4×10⁻¹² for Shirokov's
   parameters, so the bound holds with room to spare.
4. **With back-reaction the decrease probably stops.** This is a heuristic, not a derivation. Radiation lowers
   the pair's rest-frame energy E_tot below zero. By the Phase 3 1PN result d̈ = −E_tot/d², the pair then
   separates. Since A = m/d², the acceleration switches off. In the toy model the rapidity, the rest-frame
   radiated energy and the radiated 4-momentum are all finite. So M_B tends to a **finite but enormous negative
   limit** in the original frame: final rapidity ~10³ for Shirokov's numbers.

Verdict: **"M_B → −∞" holds only for the idealised eternal rigid runaway (the exact BRS solution).** For a pair
that starts from rest, a self-consistent treatment most likely gives a finite limit. That rests on a 1PN-based
toy model, not a proof.

Units G = c = 1. Scripts in this directory each write `X_output.txt`.

| script | content |
|---|---|
| `bondi_mass_check.py` | BS mass identification; Pravda–Pravdová closed-form m(u) checked against the mass-loss formula; limits |
| `linearized_news.py` | linearized retarded news of a hyperbolic ± dipole compared with the BS news |
| `radiated_momentum.py` | linearized news vs BS at 200 random points; radiated 4-momentum per unit rapidity |
| `backreaction_toy.py` | adiabatic 1PN + radiation-reaction toy model |
| `summary_numbers.py` | numbers for Shirokov's parameters; check of the u → 0 luminosity |
| `literature_notes.md` | what was read and how |

## 1. Literature: what I could and could not access

Direct PDF download from arXiv was blocked (proxy 403). I read papers only through WebFetch's model-summarised
text extraction, which is lossy: one equation came back visibly garbled. Every formula I rely on below is
**re-verified independently** in the scripts.

**Accessed (summarised text):**
- **Pravda & Pravdová, gr-qc/9806114**, "Uniformly accelerated sources in electromagnetism and gravity" (1998).
  - The BS metric is eq. (18)–(19) in cylindrical coordinates, with R = (ρ² + z² − t²)/2 and
    R_i = √((R − h_i)² + 2ρ²h_i).
  - The choice a₁ = (h₁ − h₂)²/2h₂, a₂ = −(h₁ − h₂)²/2h₁, k = 1 gives "two pairs of freely moving positive and
    negative particles" with no nodal singularity. The review gr-qc/0003067 adds: "in the case with no nodal
    singularities, negative masses occur."
  - In the small-mass limit (h₂ = h, h₁ = h + ε) the masses are m⁽¹⁾ = ε²/(2h√(2h+2ε)) and
    m⁽²⁾ = −ε²/((2h+2ε)√(2h)).
  - The news is c,u = 3ε³ u sin²θ / (2h(u² + 2h sin²θ)^{5/2}).
  - They give a closed-form m(u) (their eq. 25) and state "m as a function of u is everywhere decreasing. The
    system does radiate gravitational waves."
  - For the electromagnetic analogue they note that boost symmetry makes the rate constant in u*.
- **Bičák & Pravdová, gr-qc/9808068** (ar5iv). A second Killing vector compatible with axial symmetry and a
  smooth 𝒥 is either a supertranslation or the boost.
  - BRS news has the form K(w)/u² with w = sin θ/u.
  - Mass aspect: M = (1/2 sin θ)(w²K,w),w + L/u³.
  - Total mass: m(u) = ¼∫(w²K,w),w dθ + ½∫λ/(w²u²) dθ.
  - Null infinity is global, "though its generators are not complete."
  - No statement on whether the total radiated energy is finite.
- **Pravda & Pravdová review, gr-qc/0003067**: introduction only (quoted above).
- **Bičák & Kofroň, 1210.5870**: abstract only, nothing relevant.
- **Bičák, gr-qc/0004016** ("Selected solutions…"): outline only. BRS spacetimes are "the only known radiative
  solutions describing finite sources" and admit global smooth sections of 𝒥. Section 11 itself was not reached.

**Not accessed** (no free text): Bičák 1968 (Proc. R. Soc. A 302, 201); Bičák, Hoenselaers & Schmidt 1983
(Proc. R. Soc. A 390, 397 and 411); Bičák & Schmidt 1989 (PRD 40, 1827); Bonnor & Swaminarayan 1964.

**From memory only, not checked here.** Bičák–Schmidt showed that 𝒥⁺ of radiative BRS spacetimes admits
smooth global sections only on the half u > 0. The particle worldlines meet 𝒥 at isolated points, and the
regular part of 𝒥 is not complete. The PP mass formula is quoted for u > 0, which is consistent with this, but
I have not verified it.

## 2. The exact BS ± pair (`bondi_mass_check.py`)

**It is the Bondi dipole.**
- Particles sit at ρ = 0, z² − t² = 2h_i, so their proper acceleration is 1/√(2h_i).
- With m = ε²/(2h)^{3/2} (leading order) and proper separation d = √(2h+2ε) − √(2h) ≈ ε/√(2h), sympy gives
  **A − m/d² = 0 exactly**. This is the Bondi relation; Phase 3 found the same at 1PN.
- Mass sum: from the sympy series, m⁽¹⁾ = m(1 − ε/2h) and m⁽²⁾ = −m(1 − ε/h), with m/d = ε/2h. So
  **m⁽¹⁾ + m⁽²⁾ = +m²/d**, and the negative mass is the *lighter* one.
  - Phase 3's 1PN rigidity condition is E_tot = 0 ⇔ μ₁ + μ₂ = μ₁μ₂/d = −m²/d, with the negative mass heavier.
  - The two have the same magnitude and **opposite sign**.
  - I have not resolved this. The most likely cause is that PP's m⁽ⁱ⁾ (Curzon/Weyl mass parameters read off in
    the boosted frame) are not the 1PN harmonic rest masses μ_A. They can differ at O(m²/d) through the
    interaction or redshift factors (cf. Phase 4's E vs W weights). It could also be a sign error in my
    transcription of the extracted text. It does not affect the radiation results.
- Write a² = 2h = 1/A² and D = md = ε³/(2h)². Then the BS news becomes
  **c,u = 3 D a² u sin²θ / (u² + a² sin²θ)^{5/2}.**

**Bondi mass.**
- I took the PP closed form for m(u) with o = 3ε³/2h. Its derivative equals −½∫₀^π c,u² sin θ dθ, the
  mass-loss formula applied to their own news, to 20 digits at u = 0.1, 1, 3, 10. The (garbled) transcription
  of eq. 25 is therefore correct and internally consistent.
- Numerically m(u) = ∫_u^∞ L(u′)du′, so **m(+∞) = 0**: no residual mass.
- Limits (units D = a = 1, restored by m → D²A³·f(Au)):
  - u·m → 3/8 as u → 0⁺, and u²L(u) → 3/8, checked numerically to 6×10⁻⁷ at u = 10⁻³.
  - m = O(u⁻⁷) as u → ∞.
- The formula is odd in u: m(−u) = −m(u). Applied naively on u < 0, it would decrease from 0 at u = −∞ to −∞
  at u → 0⁻. I do not interpret the u < 0 branch, because of the 𝒥 regularity issue above.

**Where the energy goes.** For the +z pair, z = a cosh η and t = a sinh η, so emission at rapidity η reaches
direction θ at Bondi time u = a(sinh η − cos θ cosh η).
- Along θ = 0, u = −a e^{−η} → 0⁻. The particle reaches 𝒥⁺ at the corner (u = 0, θ = 0).
- The radiation from late times is beamed into θ ≲ 1/γ and arrives at |u| ~ a e^{−η}. It accumulates at u → 0.
- So the divergence m(0⁺) = +∞ is the infinite energy radiated over the infinite future of the eternal motion,
  squeezed into a finite stretch of Bondi time.
- On any regular cut u = u₀ > 0, the energy still to be radiated is finite (= m(u₀)).

**Turnaround.** At t = 0, emission reaches θ at u = −a cos θ, where c,u = −3D A² cos θ sin²θ exactly. This gives
L = ½∫c,u² sin θ dθ = **(24/35) D²A⁴ = (24/35) m⁶/d⁶**. That is exactly the Phase 3 Newtonian octupole result,
which provides an independent check.

## 3. Linearized gravity, all speeds (`linearized_news.py`, `radiated_momentum.py`)

**Source.**
- Use the signed pole–dipole tensor T^{μν} = −∂_α∫dτ S^{α(μ}u^{ν)}δ⁴(x − z(τ)), with S^{ab} = D(n^a u^b − u^a n^b)
  and zero monopole (p^μ = Σμ_A u^μ = 0 at this order).
- For any worldline in the t–z plane, S = D dt∧dz is constant. The Mathisson–Papapetrou equations with p = 0
  and constant S are satisfied identically.
- The source is therefore **conserved for hyperbolic (or any) motion**, so linearized theory is consistent. This
  is the signed source the task asked for. A single positive mass in hyperbolic motion would not be conserved
  without a strut.

**Retarded field.**
- For k = d(t − n·x), r h̄^{μν} = −4∂_u[S^{α(μ}u^{ν)}k_α/(k·u)]_ret, and c = (r/4) sin²θ h̄^{zz}.
- With the exact retarded rapidity, the single-pair news equals **−(BS news) to 1×10⁻²² at 200 random (u, θ)**:
  **c,u = −3D a² u sin²θ / (u² + a² sin²θ)^{5/2}**.
- Adding the mirror pair (z → −z) gives exactly −2 × BS.

**Unresolved: overall factor −2.** The sign is a convention. The factor 2 means either (i) PP's ε-masses or news
normalisation differ from mine, or (ii) the BS news belongs to one pair's worth of D. The single-pair linearized
luminosity matches the Phase 3 octupole independently, so I trust the linearized normalisation. Either way it
is an O(1) factor (≤ 4 in L) and does not affect any scaling below.

**Radiated 4-momentum.** Integrating the Bondi flux over the sphere at fixed emission rapidity gives
**dE/dη = (L₀/A) cosh η and dP_z/dη = (L₀/A) sinh η**, with L₀ = (24/35)D²A⁴. Checked to 12 digits at
η = 0, 0.5, 1, 2, 4. Equivalently dP_rad^μ/dτ = L₀u^μ: constant luminosity L₀ in the instantaneous rest frame,
like the Larmor formula (2/3)q²A² for a hyperbolic charge.

**Consequences for eternal hyperbolic motion released at t = 0:**
- Bondi energy E_B(t) − E₀ = −L₀t, **unbounded below**: linear in lab time t = sinh η/A, so exponential in the
  rapidity η.
- Bondi momentum −L₀(γ − 1)/A.
- The system 4-momentum P_sys = −L₀∫u dτ is past-directed timelike, with invariant mass ≈ L₀√(2t/A).

For Shirokov's parameters (m = 0.0143, d = 10, `summary_numbers.py`):

| t | η | E_B − E₀ | W₋ ≈ e^η m | E_rad/W₋ |
|---|---|---|---|---|
| 10³ | 0.14 | −5.9×10⁻¹⁵ | 0.016 | 3.6×10⁻¹³ |
| 10⁶ | 5.7 | −5.9×10⁻¹² | 4.1 | 1.4×10⁻¹² |
| 10¹² | 19.5 | −5.9×10⁻⁶ | 4.1×10⁶ | 1.4×10⁻¹² |
| 10¹⁶ | 28.7 | −0.059 | 4.1×10¹⁰ | 1.4×10⁻¹² |

Asymptotically E_rad/W₋ → L₀/(2Am) = **(12/35)(m/d)⁴**. That is below 1 whenever m/d < 1.3, i.e. throughout
the weak-field regime where the formula applies.

## 4. Back-reaction (`backreaction_toy.py`, heuristic)

**Model.** Take the instantaneous rest frame and treat the dynamics adiabatically.
- Separation: d̈ = −E_tot/d² (Phase 3, 1PN; rigid ⇔ E_tot = 0).
- Radiation: Ė_tot = −L₀(d) = −(24/35)m⁶/d⁶.
- Rapidity: η̇ = m/d².
- Scaled variables: x = d/d₀ and s = τ/T, with T = d₀³/(24m⁶/35)^{1/3} (= 5.5×10⁶ for Shirokov's numbers).

**Result.** Starting rigid (E_tot = 0):
- E_tot becomes negative, reaching −0.986 d₀³/T² = −3.2×10⁻¹¹.
- The pair separates. x grows linearly (x′ → 1.36), so A ∝ 1/τ².
- Both ∫A dτ (→ 1.93 mT/d₀² ≈ 1.5×10³ in rapidity) and ∫L₀ dτ converge.
- Starting with E_tot = ±10⁻³ (scaled) changes nothing visible. The sign is robust: radiation always pushes
  E_tot < 0, which by the Phase 3 result always means separation.

**Implication.** The radiated 4-momentum is finite, so in any fixed frame **M_B approaches a finite limit**. The
invariant loss is only ~3×10⁻¹¹, but it is boosted by e^{η_final}, which is astronomically large.

**Caveats.**
- The model uses 1PN separation dynamics at arbitrary η, justified only because it is applied in the comoving
  frame where velocities are small.
- It is adiabatic and ignores the fact that the rigid pair is unstable to mass mismatch anyway (Phase 3).
- It treats L₀ as the only loss channel.
- It is not a derivation.

## 5. Verdict and relation to Phase 4

- **Linearized and exact-BRS level (verified here):** a rigid ± pair in eternal hyperbolic motion radiates GWs
  at constant rest-frame luminosity (24/35)m⁶/d⁶. The total radiated energy is infinite, so the Bondi energy in
  a fixed frame decreases without bound, linearly in lab time.
  - In the exact BS spacetime, this infinite energy reaches 𝒥⁺ only at the corner where the particles end.
  - On the regular cuts u > 0 the Bondi mass is finite, positive and decreasing to 0.
  - So in the exact eternal solution, **no regular cut has M_B = −∞**. The divergence sits at the boundary of
    the regular part of 𝒥 (m(0⁺) = +∞ in PP's normalisation, where m(∞) = 0).
- **Phase 4 bound:** consistent, not contradicted. Loss and floor both grow like e^η. Their ratio is
  (12/35)(m/d)⁴ ≈ 10⁻¹², so radiation never comes close to the floor M_B ≥ −W₋. The floor is no more uniform in
  time than before.
- **Physical (released from rest) pair:** the luminosity is far too small to matter on simulation timescales
  (10⁻¹⁴ m per 1000 time units). The heuristic back-reaction says the runaway shuts itself off by separating,
  which gives a finite (though huge) total Bondi-energy loss.
- **Unbounded extraction by radiation is therefore not established.** It occurs only in the idealised eternal
  solution, and there only as a boundary effect on 𝒥.

**Remaining uncertainties.**
1. The factor −2 between the linearized and BS news.
2. The status of the u < 0 branch of PP's m(u) and the precise regular domain of 𝒥⁺ (Bičák–Schmidt 1989 not read).
3. The back-reaction model is 1PN plus adiabatic. A real answer needs a strong-field or 2PN radiation-reaction
   treatment, or a Shirokov-type evolution run for t ≫ 10⁶.
4. The BS mass sum m⁽¹⁾ + m⁽²⁾ = +m²/d has the opposite sign to the Phase 3 1PN rigidity condition (−m²/d).
   This is probably a mass-definition issue and remains unresolved.
5. Literature formulas were read through a lossy text extractor, and only re-verified ones are used.

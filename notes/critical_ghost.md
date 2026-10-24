# Ghost-field thresholds and Choptuik critical collapse: is it "the same coin"?

Literature check, October 2026 (branch `lit/critical-ghost`). The question is whether the knife-edge of the ghost-scalar
(Ellis) wormhole, which collapses to a black hole if perturbed one way and expands if perturbed the other way, is related
to Choptuik's threshold between dispersal and black-hole formation. The two systems are "sister equations": gravity plus a
massless scalar, with the sign of the scalar's energy flipped. Bondi, Gleiser–Dotti, the phantom stars and the mass bound
are in `notes/literature.md` and the phase RESULTS files and are not repeated here.

Follow-up: [`bh_to_wormhole.md`](bh_to_wormhole.md) works out when a black hole's horizon becomes a wormhole throat
(the sign of the energy flux through it), the Roberts family's horizons, and why a small ghost admixture fades out at
Choptuik's threshold.

Status words: **literature** (with a reference); **checked here** (by hand or sympy; scripts in the appendix);
**argued** (a derivation sketch, not a proof); **conjecture**; **unverified** (a reference detail I could not confirm).
Most papers were read through arXiv/ar5iv pages via a summarizing fetch tool, not byte-for-byte. Quoted numbers are the
papers' own.

## Summary

**The wormhole threshold is a known critical phenomenon, of type I.** The static wormholes supported by a massless ghost
scalar have exactly one unstable mode. González, Guzmán and Sarbach proved this for the whole static family [GGS09a]. It also
holds numerically in every dimension computed (n = 4–11, 20, 50, 100) [TS13], and for the charged family (proven for most
members) [GGS09c]. Shinkai and Hayward's
2002 evolutions [SH02] already measured scaling at the threshold, although they did not call it type I:
- the black-hole mass tends to a finite value, M ≈ 0.30 a (a = throat radius);
- the expanding branch tends to a Hubble rate H ≈ 1.1/a;
- the collapse time grows like 0.60 ln(1/E₀) in their null coordinate.

That last number is what type I scaling predicts from the single unstable mode. The mode's e-folding time is 0.846 a
[GGS09a, TS13], and converting to their null coordinate x⁻ = (t − l)/√2 gives 0.846/√2 = 0.598 (**checked here**; neither
paper makes this comparison). Shinkai and Torii found the same finite-mass behaviour in 5 and 6 dimensions and in
Einstein–Gauss–Bonnet gravity [ST17]. Nobody has done a careful bisection, and the published threshold masses disagree:
0.30 a [SH02] against 0.22 a [GGS09b].

**Choptuik's threshold is type II, and in the massless normal + ghost model the two live in different sectors of initial
data.**
- Choptuik's critical solution has a regular centre, is self-similar and carries no scale.
- The wormholes have two asymptotic ends and a throat, which sets their scale.
- Static solutions of any mixture of massless normal and ghost scalars cannot have a regular centre (**argued**, by Gauss's
  law, §2.4). So the regular-centre sector has no static (type I) critical solution.
- Initial data with a regular centre cannot evolve smoothly into two-ended data: a throat needs Misner–Sharp mass r/2 > 0,
  and the one exact model of wormhole formation passes through a (weak) singularity [Mae09].

A *pure* ghost field has no black-hole threshold at all. In the numerical relativity runs it always disperses
[CCP26, NRM12], and a Misner–Sharp monotonicity argument forbids trapped surfaces while the centre stays regular (§2.1,
**argued**).

**Inside the normal + ghost model, Choptuik's solution is unchanged (argued, §2.3).**
- Since T is quadratic in each field, the ghost field's coupling can be rescaled away except for its sign. So a "ghost
  fraction" is not a parameter of the equations.
- Choptuik's solution embedded in the normal field is an exact solution for either sign.
- To linear order the ghost-sector perturbations are test fields on Choptuik's background. That is the same problem as the
  imaginary-part perturbations of a complex scalar, which all decay apart from a symmetry mode [GMG96, LRR07]. So
  ghost admixture should be irrelevant at Choptuik's threshold: the same γ and Δ. No one has tested this.

**What links the two sectors in the literature, and what does not.**
- No paper found connects Choptuik's solution to the wormholes, and none fine-tunes normal + ghost data. The 2025 review
  of critical collapse does not mention ghosts or wormholes [GHMG25].
- No paper found studies critical collapse in a σ-model with a Lorentzian target. Brans–Dicke with ω < −3/2 is such a
  model, and Liebling and Choptuik explicitly left it unexplored [LC96].
- The only continuous families of exact solutions that pass from the normal to the ghost sign are two:
  - the static Fisher → Schwarzschild → anti-Fisher/Ellis family, parametrized by the signed scalar charge [BFZ11];
  - the self-similar Roberts family, which runs from black hole, through critical, dispersal and Minkowski space, to
    ghost wormhole formation ([Mae09]; sign structure **checked here**).
- Both families cross between the signs at *zero field* (vacuum), not at a threshold. Their normal-side members are not
  true critical solutions (a continuum of unstable modes [BFZ11, Fro97]).

**Verdict on "the same coin".** The sister-equation relation is exact: the ghost model is the normal model with φ → iφ. It
maps the static and the Roberts families into each other. It has no real image of Choptuik's echoing solution. The two
thresholds share a mechanism, a single unstable mode giving a codimension-one boundary with a black hole on one side. But
they are two different saddle points of one dynamical system (type II with a regular centre, type I with a throat), not one
threshold seen from two sides.

The places where a real connection could still exist, all untested:
- self-similar solutions whose fields *boost* in the Lorentzian target, the analogue of Hirschmann–Eardley's
  phase-rotating solution;
- the Brans–Dicke ω < −3/2 model;
- the outcome map of regular-centre normal + ghost data, where transient throats have been seen [NRM12].

## 1. Wormholes supported by a ghost massless scalar

### 1.1 Linear instability: one unstable mode

| system | result | ref |
|---|---|---|
| Static spherically symmetric wormholes of a massless ghost scalar, the full 2-parameter family (throat radius, product of the asymptotic masses) | **Exactly one** exponentially growing mode, proven. The operator is factorized with the zero mode 1/Ψ₀, and a nodal argument gives exactly one bound state. Time scale τ/r_throat = 0.846 (massless symmetric) → 0.590 (large asymmetry γ₁) | [GGS09a] |
| Same, in n = 4…11, 20, 50, 100 spacetime dimensions | One negative ω² in every n. ω²a² = −1.397, −2.985, −4.687, −8.290, −12.044 (n = 4, 5, 6, 8, 10); −31.6, −91.3, −191.3 (n = 20, 50, 100). Faster in higher n | [TS13] |
| Charged (Maxwell + ghost) wormholes, 4-parameter family | One growing mode, proven for sub-critical and critical members; numerical evidence for super-critical | [GGS09c] |
| Ghost scalar with a potential: Ellis–Bronnikov and an AdS wormhole | Ellis–Bronnikov: a unique bound state; the asymmetric case also has a zero-energy mode growing linearly in t (a move along the family). AdS wormhole: exactly one unstable mode | [CPS20] |
| All static V = 0 scalar-vacuum solutions, both signs | All unstable under spherical perturbations: Fisher (normal; solutions for *every* ω² < 0), anti-Fisher branches A, B and wormhole branch C. Branch A tends to Schwarzschild and its instability fades there | [BFZ11] |
| Phantom with a specific potential (M-AdS wormholes, "black universes") | Unstable, except black universes whose horizon sits at the area minimum | [BKZ12] |
| Rotating 5D Ellis wormholes | The unstable mode disappears above a critical angular momentum | [DFKKR13] |
| Neutron star with a wormhole core (ghost with quartic coupling) | Linearly unstable | [DFKK13] |

Two remarks.
- The n = 4 entries agree: 1/√1.397 = 0.846 (**checked here**).
- [GGS09a]'s proof uses the device of this repo's `phase2/dg_proof`: a static zero mode 1/Ψ₀ and its nodes count the bound
  states. The Dotti–Gleiser and GGS operators are different, but the one-mode theorem has the same structure.

*Large n (remark, my fit to three points of [TS13]):* |ω²|a² ≈ 2n − 8.6 for n = 20–100, so the growth rate rises like
√(2n)/a.

### 1.2 Nonlinear fate

- **[SH02]** (4D, dual-null code, both universes): Gaussian pulses of the ghost *or* a normal massless field. The throat
  "suffers a bifurcation of horizons." It explodes into an inflationary universe or collapses to a black hole, according
  as the total input energy E₀ is negative or positive. Normal-matter pulses (E₀ > 0) always give a black hole.
- **[GGS09b]** (4D):
  - Positive amplitude gives collapse to Schwarzschild.
  - Small negative amplitude, ε_c ∈ [−0.04, 0), gives expansion.
  - A *second* collapse window, ε_c ∈ [−0.09, −0.05], lies on the negative axis. It was scanned in steps of 0.01, with no
    bisection.
  - The expanding throat grows exponentially, with e-folding time → 0.85 as ε_c → 0, "in good agreement with perturbation
    theory". No apparent horizon forms there, and the scalar field near the throat grows.
- **[DHNS09]** (zero-mass Armendáriz-Picón wormhole, normal Φ and ghost Ψ shells):
  - A positive-energy pulse gives collapse.
  - The sign of a Ψ pulse decides the outcome.
  - The delay before collapse "becomes longer and longer" as the amplitude falls.
  - No scaling law is given.
- **[ST17]** (n = 4, 5, 6 GR and Einstein–Gauss–Bonnet, dual-null):
  - The same two fates occur.
  - The final black-hole mass "does not depend on the injected energy", which they call "a critically formed minimum mass".
  - Near threshold a temporary trapped region of constant radius appears, and the authors *guess* "first-order transition".
  - In EGB the energy threshold becomes nonzero and rises with α_GB and n.
- **Recent:**
  - [Kai26] evolves Ellis–Bronnikov and quantum-corrected Schwarzschild wormholes: same expand/collapse dichotomy, no
    fine-tuning.
  - [Shi26a] (3D, single-author preprint) gets a collapse with ℓ = 2 gravitational waves and a "phantom bounce", from
    hand-picked runs.
  - Shirokov & Nachevsky 2026 (arXiv:2610.13206, binary wormholes) is listed in `notes/novelty/search_log.md` and was not
    read.

### 1.3 Has the boundary been treated as a critical phenomenon?

Yes, by [SH02], and partly by [GGS09b] and [ST17]. Nobody did it with the standard protocol: bisection to p* in a generic
family, the lifetime against ln|p − p*|, and comparison with the linear growth rate.

**What [SH02] measured** (their eq. 3.5, Tables I–II):
- Ghost pulses: e^{x_H⁻ − 3} ∝ E₀^{−0.60}. Normal pulses: E₀^{−0.61}. Here x_H⁻ is the null coordinate where the final
  horizon forms.
- Amplitudes c_a = 10⁻¹…10⁻⁴ (ghost) and 0.5…10⁻² (normal).
- As E₀ → 0, M → 0.30 a for both kinds of pulse (M/a = 0.296 at E₀/a = 4.4×10⁻⁵).
- On the expanding side, H → 1.1/a.
- They describe this as "a new kind of critical behaviour reminiscent of that found by Choptuik".

**Reading it as type I (checked here, from their numbers).** Type I theory gives a lifetime of −(1/λ₀) ln|p − p*| + const,
where λ₀ is the single growing rate of a static critical solution [GHMG25, Sec. II D, eq. 20].
- The threshold is at E₀ = 0 for this family.
- E₀ is linear in the ghost-pulse amplitude (E₀/a = 4.4×10⁻⁵, 4.4×10⁻⁴, 4.3×10⁻³ for |c_a| = 10⁻⁴, 10⁻³, 10⁻²).
- E₀ is quadratic in the normal-pulse amplitude (6.86×10⁻⁵ → 1.71×10⁻³ for a fivefold amplitude). So in both cases
  p − p* ∝ E₀. This is consistent with the normal field's first-order perturbation decoupling from the metric (T is
  quadratic in it), and it explains why the two exponents agree.
- SH's coordinates are x^± = (t ± l)/√2 [SH02, eq. 2.26], so near the throat Δt = √2 Δx⁻.
- Their 0.60 therefore means 1/λ₀ = 0.60·√2 a = 0.85 a. The linear e-folding time is 0.846 a [GGS09a; TS13].
- The local slopes of their Table I (0.63 and 0.59 per e-fold in x⁻) bracket this. Their exponent carries two digits, so
  the agreement is at the percent level.

Two more details fit the same picture:
- the expansion rate on the other side tends to the same rate ([GGS09b]: e-folding → 0.85);
- [SH02]'s H ≈ 1.1/a is close to 1/(0.846 a) = 1.18/a.

**Is the threshold black-hole mass finite?** Yes in every study, as type I requires:
- [SH02]: M → 0.30 a;
- [GGS09b]: "the final black hole mass is universal", m_AH ≈ 0.22, with r_AH 0.51 → 0.45 as ε_c 0.07 → 0, in units of the
  initial throat radius;
- [ST17]: the final mass is independent of the injected energy in n = 4, 5, 6 (their normalization of the mass is not
  reconciled here with SH's).

The **0.30 a against 0.22 a discrepancy is unresolved.** For a single static critical solution, type I predicts one value:
the black hole at the end of the unstable manifold. Different perturbation families plus late accretion of field that has
not dispersed are a plausible cause (**conjecture**).

**Type II?** Nobody reports a power law, and none is expected: the critical solution is static and has a scale.

**A difference from textbook type I (argued).** In type I the scale usually comes from a mass scale in the field equations
[GHMG25, Sec. II D]. The massless ghost model has none. The wormhole's scale is the throat, a property of the initial data
and its topology, and the static family has a scaling zero mode. So the threshold mass is non-universal. It is a fixed
fraction of the throat radius of whichever static wormhole the data approach. Which member generic two-ended data select,
and whether [GGS09b]'s second threshold at ε_c ≈ −0.05 is the stable manifold of a *different* (rescaled or asymmetric)
member, has not been studied (**conjecture**).

## 2. Critical collapse with negative-energy or mixed matter

### 2.1 Pure ghost: no black-hole threshold

- **[CCP26]**:
  - Fully nonlinear spherical evolutions of a massless real scalar of either sign.
  - The canonical field shows Choptuik behaviour (A* ≈ 0.566 M).
  - The phantom field "always disperses". There are no trapped surfaces, naked singularities, wormholes or other stationary
    end states, even at Misner–Sharp masses near −200 M.
  - Their large-amplitude breakdown amplitude depends on the CFL number (1.43 to 6.00) and is judged numerical.
- **[NRM12]**: a phantom dilaton collapsing alone gives a non-singular spacetime at every amplitude.
- **Argued (Misner–Sharp monotonicity).** Use dual-null coordinates with ∂₊r > 0 > ∂₋r (untrapped). Then
  ∂₊E = 4πr²e^f(T₊₋∂₊r − T₊₊∂₋r) (Hayward's form, as used in [SH02]).
  - For a massless ghost, T₊₋ = 0 and T₊₊ = −(∂₊φ)² ≤ 0, so E decreases along every outgoing ray.
  - Starting from E = 0 at a regular centre, E ≤ 0 everywhere in the untrapped region.
  - A trapped or marginally trapped sphere needs 2E/r ≥ 1, and so does a throat (∇r = 0). So neither can form while the
    centre stays regular.
  - The same holds for the phantom sector of the repo's model: ε = −1 with V ≥ 0. Shirokov's sextic V is positive for
    Φ ≠ 0, since λ²/16 < μ/3. Then T₊₋ = −e^{−f}V ≤ 0 as well.
  - I did not find this argument stated in the papers read. It needs the centre to stay regular. [CCP26] and the repo's
    phase 5 runs saw no singularities.

### 2.2 Normal + ghost pairs: what has been evolved

- **No fine-tuning study of a normal + ghost massless pair was found.** The σ-model form, a target metric diag(+1, −1), was
  not found in the critical-collapse literature either.
- **[NRM12]** (Einstein–Maxwell–dilaton with a *normal* charged scalar and a phantom dilaton and/or phantom Maxwell field):
  - **Phantom Maxwell:** black holes form.
  - **Phantom dilaton:** it suppresses the singularity; for small dilaton amplitude relative to the scalar it gives a
    "dynamical wormhole".
  - **Both phantom:** non-singular, Schwarzschild-like or naked-singular outcomes; the naked singularities appear mainly at
    equal amplitudes.
  - They report transitions such as singular/non-singular near dilaton amplitude 0.2–0.25, but no critical solution or
    exponent.
- **[HY10]** (Brans–Dicke collapse of a scalar pulse, including ω < −3/2, where the BD field is a ghost in the Einstein
  frame): for ω < −1.5, oscillations of T_vv let the apparent horizon *shrink*. No critical analysis.
- **Framework for two competing massless fields** (both normal: a scalar and SU(2) Yang–Mills) [GBH19]:
  - The expected two-mode "bicritical" picture is *not* realized.
  - A quasi-DSS solution with one unstable mode connects the YM solution at large scales to Choptuik's at small scales.
  - The scalar wins at small scales because the YM perturbations of Choptuik's solution decay. §2.3 applies the same logic
    to a ghost partner.

### 2.3 Choptuik's solution inside the normal + ghost model (argued)

Model: G_ab = 8π[T_ab(φ₁) + ε T_ab(φ₂)] with □φ₁ = □φ₂ = 0, ε = ±1.

1. **The ghost "fraction" is not a parameter.**
   - φ₂ → φ₂/√|ε| maps every ε > 0 to +1 and every ε < 0 to −1, so only sign(ε) matters.
   - Equivalently, a target metric diag(1, ε) is flat for every ε. Only its signature is invariant.
   - A continuation in ε from +1 to −1 must pass through ε = 0, where φ₂ decouples from gravity. Following a fixed geometry
     through it requires φ₂ ∝ |ε|^{−1/2} → ∞.
2. **Choptuik's solution is exact for either sign.** With φ₁ = φ_Ch and φ₂ = 0, the boosts (φ₁, φ₂) → (cosh β φ₁ + sinh β φ₂,
   sinh β φ₁ + cosh β φ₂) give a non-compact family with the same metric. Along it the field amplitudes diverge as the
   direction approaches the null cone φ₁ = ±φ₂ of the target.
3. **Its linear spectrum does not see the sign.** Linearize about φ₂ = 0.
   - The φ₁ and metric sector is Choptuik's: one unstable mode, λ₀ = 1/γ, plus the gauge mode.
   - δφ₂ obeys the free wave equation on Choptuik's background, independent of ε. This is exactly the imaginary-part sector
     of a complex scalar about the real solution.
   - It contains the λ = 0 symmetry mode δφ₂ ∝ φ_Ch: the boost here, the U(1) phase for the complex field.
   - The literature says all charge-carrying perturbations of the real critical solution decay: "only one growing
     perturbation mode within the enlarged solution space" [LRR07, §3.5, summarizing GMG96].
   - In D = 5 the leading imaginary-perturbation rate is μ = −0.448 [MGM26] (as recorded in choptuik-large-d
     `matter/LITERATURE.md` §8; not re-read).
   - Identifying GMG96's charged modes at q → 0 with these free test-field modes is argued, not checked against their text.
4. **Consequence (conjecture, untested).** The black-hole threshold of regular-centre normal + ghost data near Choptuik's
   solution has Choptuik's γ and Δ. The ghost content decays or rides the boost mode. This is the same picture as for a
   complex scalar.

### 2.4 Is there a continuous family from Choptuik's solution to the static wormholes?

None was found. Five obstructions, the first four **argued**, the fifth from the literature:

- **(a) Static solutions have no regular centre.** For a constant target metric, each field obeys
  ∂_r(√−g g^{rr} ∂_rφ_i) = 0. So φ_i = Q_i u(r) + const, with the *same* harmonic function u.
  - The static map into the target is therefore a straight line, and T = (Q₁² + εQ₂²) T[u].
  - At a regular centre √−g g^{rr} → 0, which forces Q_i = 0.
  - So every static solution is a single effective scalar of sign sign(Q₁² + εQ₂²): Fisher (spacelike line), anti-Fisher
    and Ellis (timelike), or a vacuum metric with a T = 0 "stealth" field (null).
  - The regular-centre sector has no static solutions, and hence no static type I critical solution.
- **(b) Topology.** Choptuik's solution needs a regular centre and an analytic past light cone. The Ellis family has two
  ends. In regular-centre data a minimal sphere needs E = r/2 > 0 (§2.1), and [Mae09]'s exact wormhole formation passes
  through a curvature singularity at the moment of formation.
- **(c) Scale.** A static solution is not a fixed point of the scaling in τ = −ln(−t). Type I needs a scale, and here only
  the throat provides one (§1.3). The known type I ↔ type II changes happen in models with a scale in the equations:
  - Einstein–Yang–Mills: Bartnik–McKinnon against DSS, with a colored black hole at the boundary [CCB96, CHM99];
  - the massive scalar [BCG97].
- **(d) No model parameter to continue in** (2.3.1). The embedded Choptuik solution does not move with ε.
- **(e) The two exact normal ↔ ghost families pass through vacuum, not through a threshold.**
  - *Static:* Bronnikov's general V = 0 solution with ε = ±1 runs Fisher → (zero charge) Schwarzschild → anti-Fisher
    branches → Ellis wormholes. Fisher has a continuum of unstable modes. The wormholes have exactly one [BFZ11, GGS09a].
  - *Self-similar:* the metric ds² = −2du dv + R²dΩ² with R² = −uv + C₁v² + C₂u² is Roberts' family. Its Einstein tensor is
    G_uv = 0 and G_uu = v²(1 − 4C₁C₂)/(2R⁴), G_vv = u²(1 − 4C₁C₂)/(2R⁴) (**checked here**, appendix). So the scalar is
    normal for 4C₁C₂ < 1, absent (Minkowski) at 4C₁C₂ = 1, and ghost for 4C₁C₂ > 1 [Mae09].
    - At fixed C₂ > 0, C₁ runs through black hole (C₁ < 0), the critical member (C₁ = 0), dispersal
      (0 < C₁ < 1/(4C₂)), Minkowski, and ghost members (C₁ > 1/(4C₂)). The sign structure is checked here. The placement
      of the black-hole threshold at C₁ = 0 is my reading of [Mae09]'s citation of Brady 1994 for the C₁ = 0 solution;
      Brady 1994 and Roberts 1989 were not read, and the map to [Fro97]'s parameter was not checked.
    - In the ghost members R² is positive definite, and the construction forms a wormhole [Mae09].
    - The normal-side critical member has a continuous unstable spectrum, 1/2 < Re k < 1, giving β = 1 instead of the
      in-family 1/2 [Fro97]. It is not a codimension-one critical solution.
  - So in both families normal and ghost meet at zero field. The black-hole threshold and the wormhole sit at opposite ends
    of the family.
- **Change of type within a family** is known in two-field σ-models:
  - DSS ↔ CSS at κ ≈ 0.0754 through a Hopf-type bifurcation [HE97, LC96];
  - a heteroclinic/Shil'nikov-type merger at η_c ≈ 0.17 in the S³ σ-model [LTHA02, GHMG25];
  - a fold in D on the σ-model's self-similar branch (choptuik-large-d `matter/SIGMA.md` §5).

  None of these produces a static solution, because these models are scale-free too.

### 2.5 Two-field σ-models with hyperbolic, spherical or Lorentzian targets

- **Hyperbolic target** (Hirschmann–Eardley κ > 0; κ = 1 is the axion–dilaton):
  - The phase-rotating self-similar solution is critical for κ > 0.0754. Below that it has three unstable modes and
    Choptuik's solution takes over, through a Hopf-type bifurcation [HE97; GHMG25 §IV A 3].
  - Liebling–Choptuik saw the smooth DSS → CSS transition in Brans–Dicke evolutions, near λ = 2/(2ω+3) ≈ 2/3 [LC96].
  - Axion–dilaton: γ = 0.2641066 [HHS96, EHH95].
  - In D dimensions the exponent depends on the conjugacy class of the self-similar ansatz (elliptic, parabolic,
    hyperbolic) [AH20b, HK21; via choptuik-large-d `matter/LITERATURE.md` §7.3]. choptuik-large-d `matter/CSS.md` follows
    the elliptic class to D = 80.
- **Spherical target** (κ < 0): [HE97] tentatively find the self-similar solutions "ever more unstable" for κ ≲ −0.28. The
  Brans–Dicke equivalence does not apply there.
- **Lorentzian target: not studied for critical collapse.**
  - In the Einstein frame, Brans–Dicke plus a massless scalar has target metric (ω + 3/2)dσ² + ½e^{−σ}dψ², with σ = ln Φ.
    Its Gaussian curvature is −1/(4(ω + 3/2)) = −2κ for κ = 1/(4(2ω+3)) (**checked here** by hand; it matches the κ-model
    curvature quoted in choptuik-large-d `matter/CSS.md` §1).
  - For ω < −3/2 the same formula gives κ < 0. But the target is then *Lorentzian*, a constant-curvature dS₂-type surface
    with the BD field as the ghost, not HE's sphere. It is the other real form of the same complexified model.
  - [LC96] stopped at ω ≈ −1.4 and wrote: "We also anticipate interesting negative mass solutions for λ<0." [HY10]
    evolved ω < −3/2 without a critical analysis.
  - The flat normal + ghost pair is the zero-curvature member of this Lorentzian family.

## 3. Mixed positive/negative-energy configurations: thresholds and endings

Beyond `notes/literature.md` and the repo's phase 2 and phase 5 results:

- **Phantom (ε = −1) stars cannot collapse to black holes** while the centre stays regular (§2.1, argued).
  - This matches the repo's phase 5: the unstable branch migrates (δ > 0) or disperses (δ < 0).
  - It matches [Ber26]: a free ghost boson star "quickly becomes unstable, dispersing away".
  - So the phantom-star family has a stability boundary (dN/dω = 0, `phase2/radial`) but no collapse threshold.
- **Ghost matter confined in neutron stars** [Ber26] (1D NR, t ≤ 6000):
  - continuous equilibrium families;
  - unstable cases expel part of the ghost field and relax, with the star becoming more compact;
  - no black-hole collapse or threshold scan reported.
  - The predecessor treats a real ghost inside a boson star (static) [JJN23].
- **Wormhole-core stars**:
  - neutron star + ghost core: linearly unstable [DFKK13];
  - boson star with a phantom wormhole core: the throat instability persists (Dzhunushaliev et al. 2014, arXiv:1409.6978,
    already in `notes/novelty`);
  - boson stars with wormholes, static only [DMW23].
- **Black hole ↔ wormhole conversion**:
  - impulsive phantom radiation turns Schwarzschild into a static wormhole if the impulse energy satisfies |μ| < M (exact,
    [HK04]);
  - a normal scalar plus a phantom dilaton gives transient "dynamical wormholes" [NRM12];
  - the ghost Roberts solution forms a wormhole after an instantaneous weak singularity [Mae09].
- **NR with mixed-sign energy and zero total mass**:
  - a collapsing warp bubble, NEC-violating fluid [CDK24]. In the extract read, the bubble "collapses inwards and then
    dissipate[s] to infinity", with alternating positive- and negative-energy waves and an apparently stable remnant. No
    black hole was mentioned. The results section was not read.
  - The opposite-sign point-mass pair of Cao 2012 is already covered in `notes/novelty`.
- **No study fine-tunes a mixed (normal + ghost) star or cloud to a collapse threshold.**

## 4. What is open

1. **Precision type I scaling at the Ellis threshold.**
   - No bisection to p* in a generic family. The lifetime exponent is known to two digits [SH02] and has never been compared
     with 1/λ₀ in print.
   - The threshold mass is 0.30 a [SH02] or 0.22 a [GGS09b], unresolved.
   - One accurate 1+1 run along the unstable manifold of the massless symmetric wormhole would settle the mass. A bisection
     with the lifetime fitted against ln|p − p*| would test the type I law.
2. **Universality on the two-ended side.**
   - Which member of the scale (and asymmetry) family do generic near-threshold data approach?
   - Is [GGS09b]'s second threshold (ε_c ≈ −0.05) the stable manifold of another member?
   - In n > 4 and Einstein–Gauss–Bonnet, [ST17] have threshold energies and masses but no lifetime scaling.
3. **Ghost admixture at Choptuik's threshold** (§2.3). Linear analysis says it is irrelevant: same γ and Δ, ghost content
   decaying or riding the boost mode. No evolution has tested this. Unknown in D > 5: whether the decay of the ghost-sector
   (test-field) modes on Choptuik's solution persists. It is equivalent to the charge-mode question of [MGM26].
4. **The outcome map of regular-centre normal + ghost data.** Dispersal, black hole, and transient throats [NRM12] are all
   seen. Two questions are open:
   - Is there a codimension-one boundary between "black hole" and "throat/expansion" outcomes?
   - What governs it? Static solutions are excluded (§2.4a), so it would have to be dynamical or self-similar.
5. **Lorentzian-target self-similar solutions.**
   - Do regular, light-cone-analytic self-similar solutions exist whose fields boost under the homothety (the O(1,1) analogue
     of Hirschmann–Eardley's phase rotation)? How many unstable modes do they have?
   - Same questions for the dS₂-type target (Brans–Dicke ω < −3/2).
   - These are the only candidates found for a genuinely two-field critical solution mixing the signs.
6. **A continuous family from Choptuik's solution to the wormholes.** None is known. §2.4 lists the obstructions: no
   static regular-centre solutions, topology, scale, no ε-parameter. A connection, if any, would have to go through items
   4–5, and through a singular topology change.

## 5. What a two-field (normal + ghost) computation would need in choptuik-large-d

Read from `solver/dss.py`, `solver/spec.py`, `solver/matter/dss_ell.py`, `solver/matter/sigdss_solver.py`,
`matter/SIGMA.md`, `matter/SIGMA_DSS.md` and `matter/CSS.md`. Nothing was run.

### 5.1 Which equations change

The solver normalizes G_ab = 2n(∂φ∂φ − ½g(∂φ)²) with p = r(Φ + Π), m = r(Φ − Π). Since T enters linearly, two fields with
signs ε_i (ε₁ = +1, ε₂ = −1) give:
- x a_x = (a/2)[Σ_i ε_i(p_i² + m_i²) + (n − 1)(1 − a²)]. This is the only sign change in the shooting equations.
- x c_x = c(n − 1)(a² − 1). Unchanged: a massless field's contribution cancels in c = α/a for each field separately.
- p_i,τ + (x − c)p_i,x = c_x p_i − c p_i/x + nc(p_i + m_i)/(2x), and the mirror m-equation. Unchanged, one pair per field.
- The momentum constraint (a check, not imposed) and the Taylor recursions pick up the same ε_i.
  - In `centre_series`: Q2 = Σε_i(P_i·P_i + M_i·M_i).
  - In `ssh_series`: the a-coefficient aG.
- Horizon regularity per field: p_i,τ + (1 − c_x(1) − n/2)p_i = (n/2)m_i (the ℓ = 0 case of `dss_ell.py`'s condition).

### 5.2 What breaks and what does not

- **Positivity of a and the lapse survives.**
  - x(ln a)_x = ½[S + (n − 1)(1 − a²)] with S = Σε_i(p_i² + m_i²), so a > 0 while the fields are finite.
  - x(ln c)_x = (n − 1)(a² − 1) ≥ −(n − 1), so c ≥ c(x₀)(x₀/x)^{n−1} > 0, and α = ac > 0.
- **What is lost:**
  - a ≥ 1, i.e. Misner–Sharp mass ≥ 0;
  - monotonic c;
  - the Hamiltonian constraint's control of the field amplitude. Along the null directions p₁ = ±p₂ the fields do not
    gravitate at all, and the boost orbit is non-compact.
- **The self-similarity horizon gets safer.** Its transversality factor 1 − c_x(1) = 1 − (n − 1)(a(1)² − 1) *grows* where
  ghost energy dominates (a(1) < 1). The gauge mode's rate rises above 1. Nothing degenerates.
- **The gauge cannot represent a throat.** In polar-areal coordinates a minimal sphere on the slice (∇r = 0, E = r/2) shows
  up as a → ∞, exactly like an apparent horizon. With the regular-centre boundary condition, two-ended data are not in
  the function space at all. The static wormholes are therefore unreachable by this solver in principle, not just in
  practice. This applies equally to the repo's `phase5/spherical/evolve.py`, which uses the same gauge.
- **A new exact symmetry.** The boosts form a non-compact one-parameter family of solutions, so the Jacobian gets a null
  direction. A transversal condition is needed, e.g. ⟨p1₁, p1₂⟩_τ = 0. At (φ_Ch, 0) its derivative along the orbit is
  ⟨p1₁, p1₁⟩ ≠ 0.
  - If φ₂ is given the *opposite* half-period parity to φ₁ (φ₁ odd, φ₂ even), boosts leave the parity class, so the
    condition is automatic.
  - Both parity assignments are consistent truncations: φ₁ → −φ₁ and φ₂ → −φ₂ are each isometries of diag(1, −1).

### 5.3 Is a continuation in the ghost fraction set up naturally? No

- ε scales out (§2.3.1), and the embedded Choptuik solution is an exact, ε-independent point. There is nothing to continue
  in, and a Newton continuation in ε would just sit on (φ_Ch, 0).
- What *is* natural, in order of cost:
  1. **Ghost-sector spectrum of Choptuik's solution in D dimensions (cheap).**
     - Freeze the background (a, c) and solve the linear periodic problem for (p₂, m₂) alone, in both half-period parity
       sectors. This is the existing perturbation machinery (`mode_matrix` / `dss_ell.py` at L = 0) with the metric
       perturbation switched off.
     - Built-in checks:
       - the exact λ = 0 boost mode (δp₂ = p_Ch, odd sector);
       - [MGM26]'s μ = −0.448 at D = 5;
       - the D = 4 charge-mode decay [GMG96].
     - All other modes decaying in D means Choptuik's solution remains *the* critical solution of the normal + ghost model.
       This is the ghost analogue of the charge exponent, and it comes for every D.
  2. **Branch detection by continuation in D, not in ε.**
     - If a ghost-sector eigenvalue other than the boost mode crosses 0 at some D, a genuinely two-field echoing branch
       bifurcates there. That branch is what a two-field solver (5.1 + 5.2) would then follow.
     - Without such a crossing there is no branch to follow from Choptuik's solution.
  3. **Boost-rotating self-similar solutions (moderate).**
     - With φ± = φ₁ ± φ₂, one has φ₁² − φ₂² = φ₊φ₋. So T is built from ∂φ₊∂φ₋, the same structure as Hirschmann–Eardley's
       F, F̄ at κ = 0.
     - The O(1,1) model is the other real form of the complexified HE equations: φ₊, φ₋ real and independent instead of
       F̄ = F*. The CSS ansatz F = (−t)^{iω}f becomes φ± = (−t)^{±ν}f±(x), i.e. iω → ν, ω² → −ν².
     - `css_lateral.py` and `css_pert.py` already treat f and f̄ as independent analytic functions on the complex path. The
       change is the reality condition and ω → −iν.
     - The unknowns (ν, f₊(0)f₋(0)) after fixing the boost, against analyticity of f± at the horizon, give the same count as
       HE's (ω, |f(0)|). So the search is well posed.
     - For the dS₂-type (Brans–Dicke ω < −3/2) target, the same change applies to the κ ≠ 0 equations.
     - Existence is unknown (open item 5).
  4. **Anything with a throat (not natural).**
     - It needs a different gauge: proper radial distance, as in [GGS09a, TS13], or double-null, as in [SH02, ST17].
     - It needs two-ended boundary conditions.
     - Here the wormhole's linear problem in D dimensions is already solved [TS13], and it has no self-similar structure to
       share with the echoing solver.
- **On the evolution side**, the direct test of open items 3–4 is a two-field version of the repo's 1+1 code. It needs two
  real massless fields with signs ±1, bisection, and lifetime and mass fits. It stops at the first minimal or trapped
  sphere (a → ∞), so throat outcomes would need a double-null code.

## 6. Search log and limits

- **Searched (WebSearch, standard and extended):**
  - Shinkai–Hayward, GGS I/II/charged;
  - wormhole critical behaviour / type I / threshold / lifetime;
  - Torii–Shinkai and the higher-dimensional and Gauss–Bonnet evolutions;
  - Doroshkevich et al.;
  - charged and massive phantom wormholes (Bronnikov et al., Cremona et al., Batic–Dutykh);
  - phantom/ghost critical collapse; Caridi et al.; Nakonieczna et al.;
  - Brans–Dicke ω < −3/2 collapse; Liebling–Choptuik; Hirschmann–Eardley κ < 0;
  - Lorentzian-target σ-models; competing massless fields;
  - Roberts/Maeda self-similar wormholes; Hayward–Koyama;
  - mixed star + wormhole; phantom boson-star collapse; Clough et al.; recent wormhole NR (Kain, Shirokov);
  - the 2025 review [GHMG25].
- **Not covered:**
  - INSPIRE citation lists (the site needs JavaScript and the API URL was refused), so citers of [SH02] were found only
    through search and through [Kai26]'s and [Shi26a]'s reference lists.
  - [Kai26]'s bibliography (refs. 22–29: "further dynamical evolutions") could not be retrieved; none is described there
    as studying critical behaviour.
- **Not read:**
  - Shirokov & Nachevsky 2026 (2610.13206);
  - Brady 1994;
  - the GMG96 perturbation section (only its abstract and the LRR summary);
  - [MGM26] (quoted from choptuik-large-d);
  - the results sections of [CDK24] and [NRM12]'s figures.
- **Journal references** marked *(unverified)* were not checked against the publisher. All arXiv ids were opened.
- **Coverage:** results posted after about 2026-10-21.
- **AI-directed text:** none was found in the pages read.

## References

- [SH02] H. Shinkai, S. A. Hayward, "Fate of the first traversible wormhole: black-hole collapse or inflationary
  expansion", PRD 66, 044005 (2002), arXiv:gr-qc/0205041, doi:10.1103/PhysRevD.66.044005.
- [GGS09a] J. A. González, F. S. Guzmán, O. Sarbach, "Instability of wormholes supported by a ghost scalar field. I. Linear
  stability analysis", CQG 26, 015010 (2009), arXiv:0806.0608.
- [GGS09b] J. A. González, F. S. Guzmán, O. Sarbach, "… II. Nonlinear evolution", CQG 26, 015011 (2009), arXiv:0806.1370.
- [GGS09c] J. A. González, F. S. Guzmán, O. Sarbach, "On the instability of charged wormholes supported by a ghost scalar
  field", arXiv:0906.0420 (PRD 80, 024023 (2009) *(unverified)*).
- [DHNS09] A. Doroshkevich, J. Hansen, I. Novikov, A. Shatskiy, "Passage of radiation through wormholes", arXiv:0812.0702
  (IJMPD 18, 1665 (2009) *(unverified)*). A companion, arXiv:0811.2558, is withdrawn.
- [TS13] T. Torii, H. Shinkai, "Wormholes in higher dimensional space-time: exact solutions and their linear stability
  analysis", PRD 88, 064027 (2013), arXiv:1309.2058.
- [ST17] H. Shinkai, T. Torii, "Nonlinear dynamics in the Einstein–Gauss–Bonnet gravity", arXiv:1706.02070 (PRD 2017;
  volume *(unverified)*).
- [BFZ11] K. A. Bronnikov, J. C. Fabris, A. Zhidenko, "On the stability of scalar-vacuum space-times", arXiv:1109.6576
  (EPJC 71, 1791 (2011) *(unverified)*).
- [BKZ12] K. A. Bronnikov, R. A. Konoplya, A. Zhidenko, "Instabilities of wormholes and regular black holes supported by a
  phantom scalar field", arXiv:1205.2224 (PRD 86, 024028 (2012) *(unverified)*).
- [CPS20] F. Cremona, L. Pizzocchero, O. Sarbach, "Gauge-invariant spherical linear perturbations of wormholes in Einstein
  gravity minimally coupled to a self-interacting phantom scalar field", arXiv:1911.13103 (PRD 101, 104061 (2020)
  *(unverified)*).
- [DFKKR13] V. Dzhunushaliev, V. Folomeev, B. Kleihaus, J. Kunz, E. Radu, "Rotating wormholes in five dimensions",
  arXiv:1309.2448 (journal *(unverified)*).
- [DFKK13] V. Dzhunushaliev, V. Folomeev, B. Kleihaus, J. Kunz, "Mixed neutron-star-plus-wormhole systems: linear stability
  analysis", PRD 87, 104036 (2013), arXiv:1302.5217.
- [Kai26] B. Kain, "Dynamical wormholes", arXiv:2602.18231 (2026).
- [Shi26a] N. M. Shirokov, "Wormhole dynamics: nonlinear collapse and gravitational-wave emission", arXiv:2604.00071
  (2026, single-author preprint).
- [CCP26] G. Caridi, F. Corelli, P. Pani, "Challenging the weak cosmic censorship with phantom fields", arXiv:2603.13210
  (2026).
- [NRM12] A. Nakonieczna, M. Rogatko, R. Moderski, "Dynamical collapse of charged scalar field in phantom gravity",
  PRD 86, 044043 (2012), arXiv:1209.1203.
- [HY10] D. Hwang, D. Yeom, "Responses of the Brans–Dicke field due to gravitational collapses", CQG 27, 205002 (2010),
  arXiv:1002.4246.
- [GBH19] C. Gundlach, T. W. Baumgarte, D. Hilditch, "Critical phenomena in gravitational collapse with two competing
  massless matter fields", PRD 100, 104010 (2019), arXiv:1908.05971.
- [Mae09] H. Maeda, "Simple analytic model of wormhole formation", PRD 79, 024030 (2009), arXiv:0811.2962.
- [Fro97] A. V. Frolov, "Perturbations and critical behavior in the self-similar gravitational collapse of a massless
  scalar field", arXiv:gr-qc/9704040 (PRD 56, 6433 (1997) *(unverified)*).
- [HK04] S. A. Hayward, H. Koyama, "How to make a traversable wormhole from a Schwarzschild black hole",
  arXiv:gr-qc/0406080 (journal *(unverified)*).
- [Ber26] A. Bernal, V. Jaramillo, N. A. Montiel-Hernández, D. Núñez, N. Sanchis-Gual, "Gravitational confinement of ghost
  scalar fields in neutron stars", arXiv:2603.04400 (2026; PRD *(unverified)*).
- [JJN23] V. Jaramillo, Jiménez-Vázquez, D. Núñez, "Confinement of exotic matter. I. Static solutions",
  arXiv:2302.07889 (PRD 107, 064061 (2023), as recorded in `notes/novelty`).
- [DMW23] P.-B. Ding, T.-X. Ma, Y.-Q. Wang, "Study of boson stars with wormhole", arXiv:2305.19819 (JHEP 04 (2024) 033
  *(unverified)*).
- [CDK24] K. Clough, T. Dietrich, S. Khan, "What no one has seen before: gravitational waveforms from warp drive
  collapse", arXiv:2406.02466, doi:10.33232/001c.121868 (Open J. Astrophys., 2024).
- [GHMG25] C. Gundlach, D. Hilditch, J. M. Martín-García, "Critical phenomena in gravitational collapse",
  arXiv:2507.07636 (2025).
- [LRR07] C. Gundlach, J. M. Martín-García, "Critical phenomena in gravitational collapse", Living Rev. Relativ. 10, 5
  (2007), arXiv:0711.4620, doi:10.12942/lrr-2007-5.

Taken from choptuik-large-d `matter/LITERATURE.md` §11, where they were checked against arXiv/INSPIRE; their content was
not re-read here except [LC96], [HE97] and [GMG96]'s abstract:
- [GMG96] C. Gundlach, J. M. Martín-García, PRD 54, 7353 (1996), arXiv:gr-qc/9606072. δ = 0.883 ± 0.007, γ = 0.374 ± 0.001.
- [MGM26] S. Mittal, C. Gundlach, L. Martel, arXiv:2609.17207 (2026).
- [LC96] S. L. Liebling, M. W. Choptuik, PRL 77, 1424 (1996), arXiv:gr-qc/9606057.
- [HE97] E. W. Hirschmann, D. M. Eardley, PRD 56, 4696 (1997), arXiv:gr-qc/9511052.
- [HHS96] R. S. Hamadé, J. H. Horne, J. M. Stewart, CQG 13, 2241 (1996), arXiv:gr-qc/9511024.
- [EHH95] D. M. Eardley, E. W. Hirschmann, J. H. Horne, PRD 52, R5397 (1995), arXiv:gr-qc/9505041.
- [AH20b] R. Antonelli, E. Hatefi, JHEP 03 (2020) 180, arXiv:1912.06103.
- [HK21] E. Hatefi, A. Kuntz, EPJC 81 (2021), arXiv:2010.11603.
- [LTHA02] C. Lechner, J. Thornburg, S. Husa, P. C. Aichelburg, PRD 65, 081501 (2002), arXiv:gr-qc/0112008.
- [CCB96] M. W. Choptuik, T. Chmaj, P. Bizoń, PRL 77, 424 (1996), arXiv:gr-qc/9603051.
- [CHM99] M. W. Choptuik, E. W. Hirschmann, R. L. Marsa, PRD 60, 124011 (1999), arXiv:gr-qc/9903081.
- [BCG97] P. R. Brady, C. M. Chambers, S. M. C. V. Gonçalves, PRD 56, R6057 (1997), arXiv:gr-qc/9709014.

## Appendix: checks

**A1. Sign structure of the Roberts/Maeda family** (§2.4e; sympy, about 20 s):

```python
import sympy as sp
u, v, th, ph, C1, C2 = sp.symbols('u v theta phi C1 C2', real=True)
R = sp.sqrt(-u*v + C1*v**2 + C2*u**2); x = [u, v, th, ph]
g = sp.zeros(4); g[0,1] = g[1,0] = -1; g[2,2] = R**2; g[3,3] = R**2*sp.sin(th)**2; gi = g.inv()
Gam = lambda a,b,c: sum(gi[a,d]*(sp.diff(g[d,b],x[c])+sp.diff(g[d,c],x[b])-sp.diff(g[b,c],x[d])) for d in range(4))/2
G3 = [[[sp.simplify(Gam(a,b,c)) for c in range(4)] for b in range(4)] for a in range(4)]
Ric = lambda b,c: sp.simplify(sum(sp.diff(G3[a][b][c],x[a]) - sp.diff(G3[a][b][a],x[c])
      + sum(G3[a][a][d]*G3[d][b][c] - G3[a][c][d]*G3[d][b][a] for d in range(4)) for a in range(4)))
Rab = sp.Matrix(4,4,Ric); Rs = sp.simplify(sum(gi[a,b]*Rab[a,b] for a in range(4) for b in range(4)))
G = sp.simplify(Rab - Rs*g/2)
print(sp.simplify(G[0,1]), sp.factor(G[0,0]), sp.factor(G[1,1]))
```

Output: `0`, `-v**2*(4*C1*C2 - 1)/(2*(C1*v**2 + C2*u**2 - u*v)**2)`, `-u**2*(4*C1*C2 - 1)/(2*(C1*v**2 + C2*u**2 - u*v)**2)`.
A massless scalar needs G_uv = 0 and G_uu = ε(∂_uφ)²·(positive const). So ε = sign(1 − 4C₁C₂), as stated.

**A2. Arithmetic behind §1.3:**
- 1/√1.397 = 0.8461 ([TS13] n = 4) and 0.846 ([GGS09a] γ₁ = 0).
- 0.846/√2 = 0.598 against [SH02]'s 0.60 and 0.61.
- From [SH02] Table I: (8.42 − 6.98)/ln(4.38e-4/4.39e-5) = 0.626 and (6.98 − 5.64)/ln(4.32e-3/4.38e-4) = 0.585 per e-fold
  in x⁻.
- 1/(0.846 a) = 1.18/a against H ≈ 1.1/a.

**A3. Brans–Dicke target curvature** (§2.5), by hand. For A dσ² + e^{−σ}dψ²:
- K = −(1/(2√(EG)))∂_σ(G_σ/√(EG)) with E = A, G = e^{−σ}, so K = −1/(4A).
- With A = ω + 3/2 = 1/(8κ), K = −2κ.

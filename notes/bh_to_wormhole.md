# Can a black hole turn into a wormhole near Choptuik's threshold?

October 2026. A follow-up to [`critical_ghost.md`](critical_ghost.md) (the literature check, cited below as CG §n), on
two questions:
- can the black hole that forms (or fails to form) at Choptuik's threshold give rise to a wormhole;
- is it significant that a black hole seems to provide the right conditions for a wormhole?

Status words as in CG: **literature** (with a reference), **checked here** (sympy, `notes/horizon_signature.py`, log in
`notes/horizon_signature.log`), **argued** (a derivation, not a proof), **conjecture**.

## Summary

- **One equation decides black hole or wormhole** (§1, **checked here** in D = 4, 5, 6; Hayward's theorem
  [H94, H09]). A sphere that is just marginally trapped is a black-hole horizon if the null energy flowing through it
  is positive and a wormhole throat if it is negative. Both sit where the mass inside reaches the trapping value (r/2 in
  D = 4). So a black hole does provide the place for a throat. What it cannot provide is the negative energy.
- **In regular-centre collapse a throat can only appear where a black hole is forming** (§6, argued): it needs a
  marginally trapped sphere, and only positive energy can make one. A negative flux must then arrive there.
- **The exact sister family puts the wormholes on the far side of flat space from the threshold** (§2, **checked here**).
  Along Roberts' family the horizon is spacelike (black hole), shrinks to zero size and turns null at the threshold,
  vanishes (dispersal), and reappears as two timelike throats only after the field has passed through zero and become a
  ghost.
- **A little ghost field mixed into Choptuik's collapse fades out at the threshold** (§3, argued). Its energy at horizon
  formation, relative to the black hole's, falls like (p − p*)^{2|μ|γ}: exponent 0.271 in D = 4 and 0.370 in D = 5,
  from the known ghost-sector (charge-mode) rates [GMG96, MGM26]. So the first horizon is always a black-hole horizon.
  The one escape, a ghost-sector mode that grows in some D, is untested (§6).
- **Two published routes turn black holes into wormholes, and neither fits a collapse black hole** (§5): beaming
  negative energy into the eternal two-sided black hole [HK04, KH04], and Hawking evaporation, which makes the horizon
  timelike, "a traversible wormhole by any reasonable definition" in Hayward's local sense [H02]. The second is the
  closest the literature comes to "a dispersing black hole becomes a wormhole". It creates no second exit, and the
  near-threshold black holes it would act on fastest reach the Planck scale, where neither calculation holds.

## 1. One equation decides black hole or wormhole

Spherical symmetry in D = n + 2 dimensions, in dual-null form

```math
ds^2 = -2e^{-f}\,dx^+dx^- + r^2\,d\Omega_n^2 .
```

- The Einstein tensor's ++ component is G₊₊ = −(n/r)(∂₊∂₊r + ∂₊f ∂₊r) (**checked here** for n = 2, 3, 4).
- A sphere is marginally trapped where the expansion along x⁺ vanishes, ∂₊r = 0. There G₊₊ = 8πT₊₊ gives
  ∂₊∂₊r = −(8πr/n) T₊₊.
- The curve of such spheres (the trapping horizon) keeps ∂₊r = 0, so its tangent is ∂₊ + s ∂₋ with
  s = −∂₊∂₊r / ∂₋∂₊r. Its squared length is −2e^{−f}s.
- On an *outer* horizon, ∂₋∂₊r < 0 (the expansion falls going inward). Then the squared length is
  16π r e^{−f} T₊₊ / (n |∂₋∂₊r|) (**checked here**), which has the sign of T₊₊.

So the same marginally trapped sphere is
- **spacelike or null** where the null energy flux through it is positive or zero: a black-hole horizon;
- **timelike** where the flux is negative: a wormhole throat, which can be crossed in both directions.

This is S. A. Hayward's unification of black holes and wormholes (**literature**: the theorem in [H94], Theorem 2; the wormhole case in [H99, H09]). It holds in every D, and it is
the precise sense in which a black hole "provides the right conditions" for a wormhole: the black hole's horizon is the
place where a throat would be, and the only thing that decides which it is, is the sign of T₊₊ there. For a massless
scalar, T₊₊ = ε(∂₊φ)² with ε = +1 (normal) or −1 (ghost). So for the sister equations the answer is simply the sign ε:
the normal field can only make horizons, the ghost field only throats.

## 2. The exact sister family: Roberts' self-similar solutions

The one known family of exact solutions that runs from a black hole, through Choptuik-like critical behaviour, to a
wormhole is Roberts' (D = 4) [CG §2.4e; Mae09]:

```math
ds^2 = -2\,du\,dv + R^2 d\Omega^2, \qquad R^2 = -uv + C_1v^2 + C_2u^2 ,
```

with C₂ > 0 fixed. The scalar is normal for 4C₁C₂ < 1, absent (flat space) at 4C₁C₂ = 1 and a ghost for 4C₁C₂ > 1.
Its marginally trapped spheres (**checked here**):
- θ_v = 0 on the line u = 2C₁v, where R² = C₁(4C₁C₂ − 1)v². The line lies in the spacetime (R² > 0) only for C₁ < 0
  or C₁ > 1/(4C₂). Its tangent has squared length −4C₁.
- θ_u = 0 on the line v = 2C₂u, where R² = C₂(4C₁C₂ − 1)u². It exists only for 4C₁C₂ > 1, and is timelike (−4C₂).
- On the first line G_vv = −2/((4C₁C₂ − 1)v²): positive flux on the normal side, negative on the ghost side, as §1
  requires.

| C₁ | matter | marginally trapped spheres | their direction | outcome |
|---|---|---|---|---|
| below 0 | normal | one line, u = 2C₁v | spacelike | black hole |
| 0 | normal | the line u = 0, with R = 0 | null | the critical member |
| between 0 and 1/(4C₂) | normal | none | | dispersal |
| 1/(4C₂) | none | | | flat space |
| above 1/(4C₂) | ghost | two lines | both timelike | wormhole throat [Mae09] |

What the table shows:
- **The threshold is where the black hole's horizon becomes null and shrinks to zero size.** This confirms that C₁ = 0
  is the black-hole threshold, which CG §2.4e had only read off [Mae09]'s citation of Brady 1994.
- **The formula for the horizon's direction, −4C₁, does turn timelike just past the threshold, but there the sphere does
  not exist** (R² < 0). It becomes real only after the family has passed through flat space and the field has turned
  into a ghost. In the one exact family that contains both, the wormhole sits on the far side of the vacuum from the
  threshold, not next to it. Maeda's wormhole forms through a weak singularity [Mae09; CG §2.4b].
- The Roberts critical member is not a true critical solution: it has a continuum of unstable modes [Fro97; CG §2.4e].
  Choptuik's echoing solution has no known exact sister (CG: the map φ → iφ has no real image of it).

## 3. Near Choptuik's threshold, with a little ghost field mixed in

The setup of CG §2.3: a normal field φ₁ and a ghost field φ₂, regular centre, the normal field tuned near Choptuik's
threshold (parameter p, threshold p*), and a small ghost component of size ε_g. Around Choptuik's solution
(φ₁ = φ_Ch, φ₂ = 0) the ghost perturbation is a test field: it obeys the wave equation on the fixed background, and
its modes grow or decay as e^{μτ}, with τ = −ln(−t) the logarithmic scale (CG §2.3). Known values:
- μ = 0 exactly for the boost mode, δφ₂ ∝ φ_Ch;
- the leading other mode: μ = −0.448 at D = 5 [MGM26]; μ = −0.362 ± 0.012 at D = 4 ([GMG96], who write it as +0.362 with the opposite sign of τ; their
  δ = (2 + 0.362)γ = 0.883 is MGM26's relation δ_Q = (D − 2 − μ)γ).

**Argued:**
1. The solution leaves the critical regime at the scale L ∝ (p − p*)^γ, the size of the black hole that then forms
   (mass ∝ L^{D−3}).
2. By then the ghost field's non-symmetry content has shrunk by e^{μτ} at τ = −ln L, to A_g ∼ ε_g (p − p*)^{|μ|γ}.
3. On the black hole's first marginally trapped sphere, T₊₊ = (∂₊φ₁)² − (∂₊φ₂)². The ghost part is smaller than the
   normal part by a factor of order A_g², which goes to zero at the threshold. So T₊₊ > 0, and by §1 the sphere is
   spacelike: a black-hole horizon, not a throat.
4. The boost mode does not decay, but it cannot help: an exact boost (φ₁, φ₂) = (cosh β φ, sinh β φ) leaves T₊₊ unchanged
   (**checked here**). A boosted Choptuik solution forms the same black holes.
5. In energy terms: the ghost energy near the horizon, relative to the black hole's mass, falls like
   ε_g² (p − p*)^{2|μ|γ}. The exponent is 0.271 at D = 4 and 0.370 at D = 5. The closer to the threshold, the less the
   ghost field matters. Keeping it comparable would need ε_g to grow like (p − p*)^{−|μ|γ}, which leaves the
   small-admixture regime.

**So a small ghost admixture cannot open a throat at Choptuik's threshold.** The threshold keeps Choptuik's γ and Δ (CG's
conjecture), and the argument now says how fast the ghost content becomes irrelevant. It rests on μ < 0 for every
ghost-sector mode except the boost, which is known only at D = 4 and 5.

On "dispersal at the critical time": below the threshold no marginally trapped sphere ever forms, so there is no
horizon whose direction could flip. Exactly at the threshold the black hole has zero size, and Choptuik's solution
itself has no trapped spheres at all.

## 4. Shinkai and Hayward 2002 in the same language

[SH02] evolve the static ghost-field wormhole. Its throat is a double trapping horizon (both expansions zero), timelike
because the ghost flux through it is negative.
- A pulse with **positive** energy makes T₊₊ > 0 at the throat. By §1 the horizons turn spacelike and separate: the throat
  "suffers a bifurcation of horizons" and a black hole forms.
- A pulse with **negative** energy separates them the other way, into an expanding (inflationary) region.
- The boundary between the two is a type I threshold with the static wormhole as critical solution: their lifetime slope
  0.60 matches the wormhole's single unstable mode, 0.846/√2 = 0.598 (CG §1.3).

So §1 run from the wormhole side: turning a wormhole into a black hole is what ordinary (positive) energy does on its
own. The reverse, black hole into wormhole, needs negative energy, which is §5.

## 5. Other ways a black hole can become a wormhole (literature)

- **Beaming in negative energy** [HK04, KH04]. Exact constructions with impulses of negative-energy radiation:
  - The impulses are beamed in symmetrically from *both* sides of the eternal (two-sided) Schwarzschild black hole, with
    energy −M < μ < 0. When energy and timing match, the region left between them is a static traversable wormhole with
    throat radius a = 2(M + μ), smaller than the original horizon. A later positive impulse can enlarge it.
  - The two-sided topology is essential. The authors: topology change "is classically forbidden unless causal loops
    exist", and "black holes formed by gravitational collapse of supernova remnants would presumably not have the
    appropriate topology". A black hole formed from a regular centre, Choptuik's included, has one exterior.
- **Hawking evaporation** [H99, H02]. A semiclassically evaporating black hole absorbs negative energy at its horizon,
  so by §1 the horizon turns timelike. Hayward: "The evaporating object is a traversible wormhole by any reasonable
  definition" [H02]; locally these are "future wormhole horizons" [H99]. In a 2D dilaton model the end state of
  evaporation "might be a stationary wormhole" [HKL01] (an analogy).
  - This is the nearest the literature comes to the original question: a black hole that "disperses" by evaporating does
    acquire a wormhole-type horizon.
  - It is a local statement. No second exit or asymptotic region appears, and [Ter22] argues that spherically symmetric
    wormholes cannot exist in semiclassical gravity.
  - Near Choptuik's threshold the black holes are as small as one likes, so they would evaporate soonest. But they also
    reach the Planck scale, where neither the classical critical scaling nor semiclassical evaporation can be trusted. No
    paper found looks at near-critical or very small black holes this way.
- **Models that build the two ends in from the start.** A Vaidya "black bounce" turns into a wormhole when its mass drops
  below half the bounce radius, by phantom accretion or by evaporation [SMMV19]. Collapse of a fluid that violates the
  null energy condition settles into a wormhole at a nonzero minimum radius [CK21]. Neither converts a black hole formed
  by ordinary collapse.

## 6. What would change the answer

1. **Where a throat can appear (argued).** On any trapping horizon, black-hole horizon or throat, one expansion
   vanishes, so ∇r is null and the Misner–Sharp mass inside equals the trapping value (r/2 in D = 4). In data with a
   regular centre the mass starts at zero there, and only positive energy raises it (CG §2.1). So the first marginally
   trapped sphere is made by the normal field, at the place where a black hole is forming. A throat would then need a
   negative flux arriving at that sphere afterwards. This is the precise version of "a black hole provides the right
   conditions": it provides the sphere, never the negative energy.
2. **A growing ghost-sector mode.** If some ghost-sector mode other than the boost had μ > 0 in some D, the ghost content
   would grow near Choptuik's solution instead of fading, and a genuinely two-field critical solution would branch off
   there (CG §5.3). μ is known only at D = 4 (−0.362) and D = 5 (−0.448). Computing μ(D) across D on the
   choptuik-large-d background is the decisive test. It needs a test-field version of that solver (the field's two
   transport equations with the metric frozen), with the exact boost mode (μ = 0) and the D = 4 and 5 values as built-in
   checks. Not done.
3. **Large ghost content.** The arguments here are for a small admixture. Regular-centre data with comparable normal and
   ghost fields show transient throats [NRM12]; whether a sharp boundary separates black holes from throats there is
   open (CG §4.4). Item 1 says where to look: on spheres that are about to be trapped.
4. **Quantum effects** (§5).

## References

From CG (details there): [SH02], [Mae09], [Fro97], [NRM12], [MGM26] (S. Mittal, C. Gundlach, L. Martel,
arXiv:2609.17207 (2026); μ = −0.448 and Eq. (66) δ_Q = (D − 2 − μ)γ confirmed in its text for this note; its τ is
−ln(u* − u), a retarded null time).

New here (read for this note through arXiv/ar5iv pages via a summarizing fetch tool, so quotes are close but not
guaranteed byte-exact):
- [H94] S. A. Hayward, "General laws of black-hole dynamics", PRD 49, 6467 (1994), arXiv:gr-qc/9303006. Theorem 2.
- [H99] S. A. Hayward, "Dynamic wormholes", IJMPD 8, 373 (1999), arXiv:gr-qc/9805019 (journal from citations).
- [H02] S. A. Hayward, "Black holes and traversible wormholes: a synthesis", arXiv:gr-qc/0203051 (conference
  contribution; journal not found).
- [H09] S. A. Hayward, "Wormhole dynamics in spherical symmetry", PRD 79, 124001 (2009), arXiv:0903.5438. Sec. V: the
  null energy condition is violated on every wormhole mouth (T₊₊ < 0); its Eqs. (21), (24) give the slope formula of §1.
- [HK04] S. A. Hayward, H. Koyama, "How to make a traversable wormhole from a Schwarzschild black hole", PRD 70, 101502
  (2004), arXiv:gr-qc/0406080 (journal from [H09]'s reference list).
- [KH04] H. Koyama, S. A. Hayward, "Construction and enlargement of traversable wormholes from Schwarzschild black
  holes", PRD 70, 084001 (2004), arXiv:gr-qc/0406113, doi:10.1103/PhysRevD.70.084001.
- [HKL01] S. A. Hayward, S.-W. Kim, H. Lee, arXiv:gr-qc/0110080 (2D dilaton gravity; journal *(unverified)*).
- [GMG96] C. Gundlach, J. M. Martín-García, "Charge scaling and universality in critical collapse", PRD 54, 7353 (1996),
  arXiv:gr-qc/9606072, doi:10.1103/PhysRevD.54.7353. Sec. VI: λ₁ = −2.674 ± 0.009, μ₁ = 0.362 ± 0.012 with
  τ = ln(t/r₀), γ = −1/λ₁, δ = −(μ₁ + 2)/λ₁.
- [SMMV19] A. Simpson, P. Martín-Moruno, M. Visser, CQG 36, 145007 (2019), arXiv:1902.04232.
- [CK21] S. Chakrabarti, S. Kar, PRD 104, 024071 (2021), arXiv:2106.14761.
- [Ter22] D. R. Terno, arXiv:2203.03770 (journal not found).

# Can a black hole turn into a wormhole near Choptuik's threshold?

October 2026. A follow-up to [`critical_ghost.md`](critical_ghost.md) (the literature check, cited below as CG §n), on
two questions:
- can the black hole that forms (or fails to form) at Choptuik's threshold give rise to a wormhole;
- is it significant that a black hole seems to provide the right conditions for a wormhole?

Status words as in CG: **literature** (with a reference), **checked here** (sympy, `notes/horizon_signature.py`, log in
`notes/horizon_signature.log`), **argued** (a derivation, not a proof), **conjecture**.

## Summary

PLACEHOLDER

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

This is S. A. Hayward's unification of black holes and wormholes (**literature**, REFS). It holds in every D, and it is
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
- the leading other mode: μ = −0.448 at D = 5 [MGM26]; at D = 4, MU4.

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
   ε_g² (p − p*)^{2|μ|γ}. The exponent is EXP4 at D = 4 and 0.370 at D = 5. The closer to the threshold, the less the
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

## 5. Other ways a black hole can become a wormhole

PLACEHOLDER

## 6. What would change the answer

PLACEHOLDER

## References

PLACEHOLDER

# Which "universe from collapse" claims can a computation here refute? (#211)

November 2026. The papers are recorded in [literature.md](literature.md), section "Universe from collapse"; what was
read of each, and the verification of each reference, is in negative-mass-raw `I211-rpi5/README.md`. Context:
`critical_ghost.md` §1.1–1.3, `black_universe/README.md` (#212), `wormhole_evolution/README.md` (#39). Status words as
in AGENTS.md. Nothing here says that a universe forms; each item is a claim that a stated computation can fail.

## Result

- Of the proposals read, two have a classical solution that this project already evolves or counts: the ghost wormhole's
  expanding branch ([SH02], [GGS09b]) and the phantom black universes ([BF06], [BKZ12]). Their open claims are the end
  state of the expansion (#213), the nonlinear fate of the one-mode black universes (#217) and the non-radial and M-AdS
  counts (#218).
- The central claim of [BF06], [BMD07] and [BBS12], a universe from phantom-dominated collapse, has no formation
  calculation in any paper read, and no issue covers it. Three tasks below are covered by no issue (T1, T2, T4), and T3
  extends #42.
- The proposals of Smolin, Popławski and the modified-gravity collapses rest on physics outside these tools. For the
  limiting-curvature and S-brane gluings only the junction algebra is checkable (T4).

## 1. The proposals

Details and quotations are in `literature.md`. "Gives" is what the paper contains: S a solution, F its formation by
collapse, Stab its stability.

| proposal | claim | matter and energy condition | gives | tools here |
|---|---|---|---|---|
| [SH02], [GGS09b] | negative-energy pulse: throat "explodes to form an inflationary universe"; "expands forever" | massless ghost scalar; NEC fails (R_kk = −2(k·∇φ)² ≤ 0) | S, numerical evolution ([SH02] to x⁻ = 10); no analysis of the end state | #39 code, #213 |
| [BF06], [BMD07], [BBS12] | "our Universe could be created from a phantom-dominated collapse in another universe" | phantom scalar with V(φ), NEC fails where r″ > 0; BBS12 adds a magnetic field | S only; the collapse claim is stated as speculation | #212 code, evolution to be built |
| [BKZ12] | all unstable radially except c = −1 | the [BF06] family | Stab (linear, radial and axial) | done in #212 |
| [SV19], [BW22] | a "bounce into a future incarnation of the universe" | a metric; source found by [BW22]: phantom scalar plus magnetic nonlinear electrodynamics; NEC fails everywhere except horizons | S; no F, no Stab | #212 method, extended |
| [PPS25], [CK21] | collapse through a bounce, thin shell or imperfect fluid | no matter model; the shell or fluid is read off the geometry | F by a matching or an ansatz, not from a field theory | none |
| [FMM89], [FMM90], [EB01], [BHR21] | interior glued to de Sitter at a spacelike surface; a new universe | limiting curvature; S-brane with tension, "violates the usual energy conditions" | matching only; no F, no Stab | junction algebra (T4) |
| [Smo92], [Smo06] | "Black hole singularities bounce and evolve to initial states of expanding universes" | none stated | none; an ensemble hypothesis with observational predictions | none |
| [Pop10], [Pop16] | each black hole contains a new universe | Einstein–Cartan spin–torsion repulsion | S in the cosmological stage; no F computed | none |
| [MSOF23], [BCHM25] | collapse gives an expanding baby universe or a bounce into a new universe | Palatini f(R) boson stars; higher-curvature gravity with dust | F, numerical or analytic, in other theories | none |

## 2. Ranking

Order: bearing on the universe claim (formation above end state above stability above existence), whether it can fail,
cost with the tools here, and whether an issue already covers it. Cost: S under a day, M days, L new code.

| rank | claim a computation could refute | from | bearing | cost | covered by |
|---|---|---|---|---|---|
| 1 | the expanding branch tends to a de Sitter-like state: H constant, w → −1, an anti-trapped region, masses conserved at the ends | [SH02], [GGS09b] | high | M | #213 (open) |
| 2 | a universe forms from phantom-dominated collapse: a trapped surface forms from regular-centre data, with a bounce inside | [BF06], [BMD07], [BBS12] | high | L | none: T1 |
| 3 | a black universe with one growing mode is a critical solution, and its two outcomes say which side expands | implied by [BKZ12] and #212 | medium | M to L | #217 (open) |
| 4 | the black bounce and the magnetic black universes are radially stable or not | [SV19], [BW22], [BBS12] | medium | M | none: T2 |
| 5 | non-radial and M-AdS modes of the black universes | [BKZ12] (not computed there) | medium | M | #218 (open) |
| 6 | SH02's sign rule (E₀ < 0 expands), its M ≈ 0.30a and H ≈ 1.1/a, and the second window of [GGS09b] | [SH02], [GGS09b] | medium | M | #40, #41, #42, #150 (open); T3 extends #42 |
| 7 | the junction algebra: positive tension is needed for expansion, and the shell violates the NEC | [BHR21], [FMM90] | low | S | none: T4 |
| – | the family of [BF06] solves its equations | [BF06] | – | – | done: #212 S1, S2 |
| – | the exact black universe is asymptotically isotropic | [BF06], [BKZ12] | – | – | proved (by hand), below |
| – | the unstable-mode rates and the one stable member | [BKZ12] | – | – | done: #212, 20 of 20 agree |

Two items hold already and cannot fail within the family.
- Printed equations. [BF06] (8) halves A′(h) (by hand from its (7), 40 digits in #212), and [BKZ12] (36)–(38) disagree
  with the derivation while (40), (41) agree (#212). Both are in `literature.md`.
- Isotropy. From [BF06] (10), (13): B = A/r² → B₋ = (2bc − πρ₀)/(2b³) as ρ → −∞, and with 2bc = −πρ₀ this is −πρ₀/b³ < 0
  for ρ₀ > 0, so the far region is a T region. With r ≈ |ρ| both scale factors, r and √|A|, grow like |ρ|, and a/r =
  √|B| tends to √|B₋| with a deviation O(ρ⁻²) from the 1/r² term of (10). The shear to expansion ratio therefore tends
  to zero, and the expansion rate tends to H = √(−B₋) = √(2|c|)/b (*proved* by hand, for the exact background only). It
  does not show that perturbed or collapse-formed interiors isotropize, which [BF06] attributes to "particle creation".

## 3. Remarks on #213 and #217

These two issues cover rank 1 and rank 3. Three points from the reading, for the coordinator to take or leave.

- **What [SH02] measure as H.** Their fit is r/a = 1 + b₄ exp(H(τ − b₅)) on the diagonal (§III B), so H is the exponent
  of the departure r/a − 1, the same form as the single growing mode. For small E₀ they find Ha → 1.12 (Table I); the
  linear rate is Wa = 1/0.846045 = 1.182 at D = 4 (`wormhole_largeD/w2_output.txt`, via #39). A match to 5 % in the
  small-E₀ limit shows the fit follows the linear mode there, and does not by itself show a de Sitter state. A test that
  can fail: ln(r/a) against the throat observer's proper time keeps its slope once r/a ≫ 1 (r on the diagonal reaches 25
  a in #39's probe), or it does not. #39's probe has H on the diagonal rising 0.75 → 2.9 and the diagonal's proper time
  saturating; neither is interpreted.
- **What w can be.** For a homogeneous massless ghost field, T_ab = −(∂_aφ ∂_bφ − ½ g_ab (∂φ)²) gives ρ = p = −½φ̇², so
  w = +1 and ρ < 0 (*proved* by hand for the homogeneous case). A w → −1 would then need the gradient terms of an
  inhomogeneous region, or the reading of the end state as de Sitter fails. This is a reason for #213 to fix its
  observer and its definition of w before reading H(t), as #39's README already says; the statement is a *conjecture*
  for the inhomogeneous region.
- **Update from #213:**
  - The slope test is done and not met: d ln r/dτ crosses W near r/a = 6.5 and reaches 2.75 W at r/a = 50 (D = 4, tuned
    member). #39's probe values (H = 0.75 to 2.9, proper time saturating) are the same behaviour, with 1/H linear in τ.
  - The conjecture above for the inhomogeneous region is settled for the mirror-symmetric throat. With a purely radial
    gradient the isotropic average is w̄ = (1 − n)/(1 + n), −1/3 at D = 4 (*proved*). So a de Sitter reading fails there.
  - The masses at the two ends change by the ghost flux through them and are not conserved, and the anti-trapped region
    reaches the grid's edge. That bookkeeping supports no statement about a universe (`wormhole_evolution/expansion.md`).
- **Where the universe would be.** The expanding region of #213 has both ends asymptotically flat. [SH02] call the
  result a universe and say the two trapping horizons "become the cosmological horizons"; the mass bookkeeping in #213's
  acceptance is the test of how that region relates to the ends. #217's evolution covers the interior (T) region too,
  which #212 did not.

## 4. Proposed tasks

Each can fail. Status words as above; none of them may say that a universe forms.

### T1. Does phantom-dominated collapse with a black-universe potential form a trapped surface, and a bounce inside?

- **What.** Evolve regular-centre data, a Gaussian phantom shell of amplitude p, with the potential V(φ) of [BF06] (14)
  for c = −1.5 and c = −0.3 (b = 1, the members of #217). Scan p. Classify each run as dispersal, trapped surface
  formed, or neither, with the horizon and anti-trapped-surface finders of #39's code. For the runs that trap, continue
  inside: does r reach a positive minimum (a bounce) or go to zero?
- **Why.** [BF06], [BMD07] and [BBS12] propose our universe from phantom-dominated collapse, as "speculate" ([BF06]) or
  "make possible" ([BBS12]). No paper read computes it. `critical_ghost.md` §2.1 argues that a phantom scalar with V ≥ 0
  and a regular centre cannot trap a surface; [BF06]'s no-go (B) says a regular-centre solution with m > 0 needs a
  potential of alternating sign, so whether V(φ) ≥ 0 on the range of φ that the data visit is not settled here. The
  result either supports or fails the formation step of the proposal for this family.
- **Acceptance.**
  - Planted: with V = 0 the code never traps a surface at any p (the result of [CCP26] and the §2.1 argument); a trapped
    surface there means the code or the argument is wrong.
  - Control: ε = +1, V = 0 gives a threshold and Choptuik's γ = 0.374 at D = 4 within a stated fit error (the control of
    #45).
  - Planted: the static member of #217 stays static to the truncation error over a declared window with the potential
    term.
  - The outcome map at two resolutions that agree, and each trapped-surface time and each minimum of r reported with its
    definition. Each check has a `--mutate=` that makes it fail, computed in a development run first.
- **Depends.** #39's code with the potential term of #217 and a regular centre (not built; #45 and #46 plan a centre for
  normal plus ghost fields). Related: #45, #46, #217.
- **Covered by.** None.

### T2. Radial stability of the black bounce with its source, and of the magnetic black universes

- **What.** Linearize about the static solution as #212 did, with a phantom scalar with V(φ) plus a magnetic field: the
  nonlinear electrodynamics L(F) of [BW22] for the Simpson–Visser family (a/m from the wormhole through the one-way
  wormhole a = 2m to the black bounce, at least 10 values), and a Maxwell field for the solutions of [BBS12] (at least
  10 charges). Count the growing radial modes with #212's shooting and its spectral check.
- **Why.** [SV19] describe a "bounce into a future incarnation of the universe" and prove nothing about stability.
  [BW22] say stability "cannot be formulated until we know their precise dynamic context" and that it may depend on the
  source; this source is now known. [BBS12] ask whether electric or magnetic fields stabilize the unstable members. In
  spherical symmetry a magnetic field carries no new radial degree of freedom (its flux is fixed by the area), so a
  Schrödinger problem for the scalar should remain (*conjecture*; to be checked in the derivation).
- **Acceptance.**
  - Planted: Simpson–Visser at m = 0 is the Ellis wormhole, so the code returns W²a² = 1.39705243433 (D = 4,
    `w2_output.txt`) to #212's tolerances; and zero charge returns #212's table.
  - Counts at two resolutions that agree, with the method and the range searched, and an independent spectral count as
    #212's S6.
  - The linearized equations checked symbolically against the Einstein and field equations (residual 0), as #212's S4.
  - Each check has a control that changes the count, computed before it is declared.
  - The result is stated for this source only ([BW22]: other sources of the same metric may differ).
- **Depends.** #212 (code, method). Related: #218.
- **Covered by.** None.

### T3. Does SH02's sign rule hold at large amplitude?

- **What.** In #39's code, scan a one-sided ghost pulse over a wide range of amplitude of both signs at fixed width, as
  [SH02] did for |c_a| ≤ 0.1 (c_b ≥ 3). For each run record E₀ with [SH02]'s definition (Bondi energy at x⁺ = 20 minus
  the static value) and the outcome. Look for a collapse with E₀ < 0, as in the second window of [GGS09b], ε_c ∈ [−0.09,
  −0.05].
- **Why.** [SH02]: "positive E₀ will cause collapse to a black hole, while negative E₀ will cause explosion". [GGS09b]
  find a second collapse window on the negative-amplitude side, with data that leave the ADM masses unchanged (their
  §II), so the two energies are not comparable as defined. #41 reconciles the mass definitions, and #42 asks about the
  second window; this isolates whether the sign rule itself is a small-amplitude statement.
- **Acceptance.**
  - Outcome against E₀ at two resolutions that agree, over a stated amplitude range.
  - Either the rule holds over that range, or a run with E₀ < 0 that collapses is reported with its E₀.
  - Control: normal-field pulses always collapse ([SH02] Table II).
  - A `--mutate=` that flips the sign of E₀ makes the check fail.
- **Depends.** #39 (done). Extends #42; uses #41's definitions.
- **Covered by.** Partly #42.

### T4. Junction algebra of the S-brane and limiting-curvature proposals

- **What.** In sympy, derive the Israel conditions on a spacelike hypersurface r = r₀ inside Schwarzschild to (i) a
  homogeneous anisotropic cosmology and (ii) de Sitter. Compare with [BHR21] (4)–(5), which print the matching and say
  that "a positive S-brane tension is required in order to get expansion". Compute the surface stress and its NEC and
  WEC status for each sign of the tension.
- **Why.** [BHR21] say the S-brane "violates the usual energy conditions"; [FMM90] attach Schwarzschild to de Sitter "at
  some spacelike junction surface" with Israel's method, and the abstract does not say what the layer's stress is. #212
  found printed equations that do not match a derivation; the same check is cheap here.
- **Acceptance.**
  - The printed [BHR21] (4)–(5) reproduced, or the discrepancy reported with the derived forms.
  - The sign of the tension that gives expansion, and the NEC sign of the shell, each by hand and symbolically.
  - A mutation that flips the sign of the jump of K fails the check.
  - [FMM90]'s text is not available here (publisher HTTP 403); the part on it needs the coordinator to supply the paper.
- **Depends.** None.
- **Covered by.** None.

## 5. What these tools cannot test

- [Smo92], [Smo06]: cosmological natural selection is a hypothesis about an ensemble of universes with predictions about
  particle-physics parameters and astrophysics. It names no bounce mechanism ([Smo06] cites Wheeler and a 2005 review by
  Bojowald). A classical bounce such as a black universe would be an example of the hypothesis, not a test of it.
- [Pop10], [Pop16]: Einstein–Cartan gravity with a spin fluid is a different theory from the one coded here.
- [MSOF23], [BCHM25]: Palatini f(R) and higher-curvature collapses are in other theories. [MSOF23] read the baby
  universe in the f(R)-frame metric, while the Einstein frame gives an ordinary black hole.
- [SH02]'s "mechanism for inflating, to macroscopic size, a Planck-sized wormhole formed in space-time foam": the foam
  is outside classical evolution. The classical part is rank 1 and rank 6.
- The observational bound on a Kantowski–Sachs universe that [BF06] cite, and the Planck-scale or string-scale inputs
  (the limiting curvature, the S-brane tension).

## 6. Search limits and flags

- The papers beyond the issue's list came from [BHR21]'s reference [3] and from INSPIRE citer lists of [BF06] and [SV19]
  (title words stability, collapse, formation, dynamics; raw `inspire_citer_searches.txt`). No wider search was made,
  and the citers of [BKZ12], [BW22] and [GGS09b] were not listed. Titles only for Lobo et al. (arXiv:2003.09419) and
  Franzin et al. (arXiv:2201.01650, 2207.08864, 2310.11990), which bear on dynamic black bounces and on perturbations of
  regular black holes and wormholes; T2 should read them first.
- Not read: [FMM89], [FMM90] and [Smo92] full texts (abstracts only); [BKZ12] §III and its quasinormal-mode computation;
  [Smo06] §4; [BHR21]'s companion arXiv:2103.02842.
- Bibliographic discrepancies, none changing an entry: [BW22]'s arXiv page prints "Phys. Rev. D 15" (INSPIRE: 105);
  [Hay06]'s arXiv title says "non-singular", the journal's "regular"; [Smo92] ends on page 192 on INSPIRE, 191 in
  [Smo06].
- `critical_ghost.md`'s summary calls [SH02]'s H ≈ 1.1/a "a Hubble rate" of the expanding branch, and its §1.3 notes
  that it is close to 1/(0.846 a). Section 3 above adds that [SH02]'s H is the exponent of the departure of r/a on the
  diagonal (1.12 against the linear 1.182). Not edited here.

## Coordinator's review

Under review: `wip/i211-universe-lit` at eac1f17, raw `I211-rpi5` at 77cb32e. Review folder
`negative-mass-raw/I211-mac-review`, whose README holds the criteria, committed before the first run of `ind211.py`.

- **References.** A fixed sample of seven of the section's 25 arXiv ids was checked against arXiv's API: every entry
  carries arXiv's title and first author (ind211 K1, run 2). Run 1 failed on the review's own handling of "González"
  (amendment A1).
- **[SH02]'s H.** Read from [SH02]'s TeX source: the radius "can be expressed with a function r/a = 1 + b₄ exp(H(τ − b₅))",
  with H = 1.19, 1.13 and 1.12 in its table as the amplitude falls. The note's reading holds, and `critical_ghost.md`
  §1.2's "Hubble rate" now carries an update saying so (ind211 K2).
- *Added by the coordinator:* the claim that a homogeneous massless ghost has w = +1 is elementary, ρ = p = εφ̇²/2 for
  either sign, and stands as *proved*.
- **The ranked claims (§2) and the tasks (§4).** T1, T2 and T4 are filed as #222, #223 and #224. T3 (the sign rule at
  large amplitude) is added to #42, which it extends. #213, #217 and #218 already cover claims 1, 3 and 5.
- **[FMM89] and [FMM90].** These are read from their abstracts only (publisher HTTP 403). The review could not reach
  the texts either; their entries say so.

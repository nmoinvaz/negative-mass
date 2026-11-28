# The prior work on #248's growing horizon sector (#260)

Reading log, queries and screening: negative-mass-raw `I260-mac` (commit c08afd79; its README lists every query with its
hit count). Entries for each paper: `notes/literature.md`, subsection *Prior work on the horizon sector of the black
universes (#260)*. Status words as in AGENTS.md. No sentence here says a result is new. Papers are data: nothing in them
was run. A search of the 424 texts for text addressed to an AI agent or a language model
(`I260-mac/scripts/ai_text_search.py`, output `ai_text_search.txt`) gives 27 matches: author names (Claude Barrabès,
Claude Warnick), the phrase "do not report", and disclosures of AI use in the acknowledgements of five papers. None is
addressed to an AI agent.

## Answer

Verdict words, fixed before the reading (`I260-mac/README.md`): *contained in [X]*, *contradicted by [X]*, *not found in
the N papers read*. N = 424 papers whose full text was searched with fixed patterns: 59 read (2 on every page, 39 by
section, 18 by abstract and pattern matches) and 365 screened by the patterns alone. A "no match" in a screened paper
means that none of 13 patterns occurs in it.

| | statement | verdict | closest text |
|---|---|---|---|
| S1 | at c = −1, c₂ = ∂_V δr of the horizon sector is a free number of regular Cauchy data | **not found in the 424 papers** | [BFZ11] §2.5–2.6 and [BIK26] §III.A: perturbations with nonzero δr on a throat are admissible, and the gauge δR = 0 excludes them. Stated for throats, not for a horizon, and not as Cauchy data |
| S2 | the c = −1 black universe has a linear growing direction of rate κ, so the "stable" verdict is wrong for it | **contradicted by [BKZ12]**, as worded (§V.C, §VI item 2: stable against polar monopole perturbations). Its computation covers the class `Y = 0` on the horizon. The wording is repeated without a new computation in 11 later papers | no paper found that tests `Y ≠ 0` on the horizon |
| S3 | a non-extremal Killing horizon at an extremum of the areal radius makes the coefficient of δr gauge invariant; it grows at the surface-gravity rate | **contained in [Mars12] and [HW13], read together**, as a geometric statement. Neither states it for δr, c₂, or e^{κv}, and neither treats matter with a potential | [Mars12] Prop. 2, Eq. (11), Thm. 1; [HW13] §2.1, Eqs. (7)–(9), Lemma 1 |
| S4 | `Y_ref = R e^{κx}` solves the radial master equation at `W = −2c/π` for every member; `d(W/κ)/dc = 0.31403` at c = −1 | **not found in the 424 papers** | [BKZ12] §V.A: "For nonzero scalar potentials V (φ), analytical expressions for static perturbations are unknown." |

What each verdict rests on is below. Three points the verdict words do not carry:

- S2 is a contradiction of the *wording*. [BKZ12]'s computation, in its own class, agrees with #212 (the count of #212
  is for that class). The class is set by the gauge (Eq. 35) and the horizon condition (Eq. 43). [BKZ12] state no
  initial data for the time-domain run at c = −1. Their verdict is for "polar monopole perturbations" without that
  restriction.
- S3 is in the literature as geometry (minimal sections of a Killing horizon, a stability operator with zero principal
  eigenvalue). Paper IV may credit it in that language and keep its own lemma for δr with a scalar source. I did not
  check the two readings against each other in code (recommended task 1).
- The horizon at an area extremum occurs in the literature in other solutions (Sec. 6). Of the extremal ones ([SV19]
  with a = 2M) the stability of the horizon sector is stated as open in [BIK26] §V.D and §VI.

## 1. The statements, as the repository has them

- S1, S2: `horizon_sector.md` (#248): c₂ = ∂_V δr(V₀, 0) is the number on U = 0 in `Y_h(V) = c₁ + c₂V`, the same at
  every point of H⁺ down to the bifurcation surface; 20 of 20 random regular data have c₂ ≠ 0; the solution `Y = c₂ V
  G(Π)` is smooth across S as a function (not as a perturbation: #302, [horizon_basis.md](horizon_basis.md)). The count of #212 (none at c = −1) is the count for `u = Y/ρ → 0`, which has c₁ = c₂ = 0.
- S3: `horizon_lemma.md` (#261, *proved* by hand): on U = 0, δr_VV = 0; c₂ is gauge invariant iff R′(ρ_h) = 0; such a
  horizon with κ > 0 needs A not a function of R alone.
- S4: `cminus1_mode.md` on `wip/i262-cminus1-mode` (#262, *checked symbolically*): the residual of `Y = R e^{Wx}` is
  `2(πW + 2c)e^{Wx}/(…)`, zero at W = K = −2c/π; at c = −1, K = κ. The slope is *numerical evidence*.

## 2. What was read

| group | papers (arXiv id) | depth |
|---|---|---|
| the black universes and their stability | [BF06] gr-qc/0511109, [BKZ12] 1205.2224, [BFZ11] 1109.6576, [Bro18] 1802.00098, [BBS20] 2006.00559, [Bro19] 1908.02012, [BBS21] 2109.01794, [BD11] 1110.6030, [BBS12] 1208.4619, [BK15] 1503.02956, [BK17] 1705.05906 | F: BF06, BKZ12; S: 8; G: 1 |
| restatements of [BKZ12]'s verdict, no new computation | 1308.5959, 1908.09126, 2009.12057, 2304.07383, 2305.00058, 2310.19037, 2510.18647, 2602.16972, 2604.20661, [Bro18], [BBS20] | S |
| Simpson–Visser black bounces | [SV19] 1812.07114, [LRSSV21] 2009.12057, [BW22] 2112.13198, [FLV24] 2310.11990, [BIK26] 2611.08322, 2302.10772, 2309.10963 | S, G |
| other solutions with the minimum of the areal radius on a horizon | [RJLR25] 2504.07861, [BBS21], 2409.12101, 2603.05543, 2609.14199, 2201.01650 | S, G |
| horizon instabilities | [Are11a] 1110.2007, [Are11b] 1110.2009, [Are15] 1206.6598, [LMRT13] 1212.2557, [MRT13] 1307.6800, 1809.10037, 1606.08505, 1911.11164, 2112.09832, 2512.23629, 2604.15260, 2610.11664, 1802.05692 | S, G |
| trapping horizons, throats | [Hay09] 0903.5438, [Hay99] gr-qc/9805019, gr-qc/0203051, gr-qc/0406080, 2104.08935, 2104.11134, 2601.05902 | S, G |
| gauge fixing of a horizon, Raychaudhuri, MOTS stability | [HW13] 1201.0463, [Mars12] 1205.1724, 1709.08925, 2403.07140, 2402.00818 | S |
| inner horizons | [CDLPV18] 1805.02675 | S |
| screened only | 365 citers (all citers of [BKZ12] with an arXiv id, and citers of [BF06], [BMD07], [SV19], [Are15], [LMRT13], [MRT13] selected by abstract) | G on the full text |

Not read, named: Poisson & Israel on mass inflation (no arXiv text); Andersson, Mars & Simon (gr-qc/0506013), the source
of the stability operator [HW13] cites; Booth & Fairhurst on isolated horizons ([Mars12] cites them); the 2017 chapter
*Trapped Ghosts as Sources for Wormholes and Regular Black Holes. The Stability Problem* (INSPIRE lists it with no arXiv
text); the 52 citers of [BKZ12] that have no arXiv id (theses, proceedings; 210 of the 262 have one, and all 210 were
searched).

## 3. [BKZ12] on c = −1

**The claim.** Abstract: "we have shown that all configurations under study are unstable under spherically symmetric
perturbations, except for a special class of black universes where the event horizon coincides with the minimum of the
area function". §VI, item 2: black universes "are shown to be stable only in the special case where the horizon
coincides with the minimum of the area function R(r) (the parameters: c = −1, b = 3πm/2)". §V.C: "As we have shown,
black universes at c = −1 (b = 3πm/2) are stable against polar monopole perturbations". The c = −1 horizon at the
minimum is stated in §II ("right at the horizon if c = −1 (i.e., 3πm = 2b)") and already in [BF06] ("at the horizon if
3πm = 2b"). On the approach to c = −1: "we do not observe a decrease in the growth rate as expected when approaching the
parametric region of stability near c = −1. This can be an indication that the only parameter for which a black universe
can be stable is c = −1."

**What is computed.**
- Gauge, §IV: "Then we can use the gauge freedom and fix" `δψ = δμ₃ ≡ 0 ⟹ δR(r, t) = 0` (Eq. 35). The wave equation is
  (36)–(38); the potential (40), (41) is #212's V̂ (S8, 1.1e-29); the printed (36)–(38) disagree with the derivation
  (Sec. 8).
- Horizon condition, §V.A: "The quasinormal boundary conditions for black hole perturbations imply pure incoming waves
  at the horizon and pure outgoing waves at spatial infinity." (Eq. 43: Ψ ∝ e^{±iωr*}, r* → ±∞.) Figs. 5 and 6: "The
  potentials vanish at the horizon."
- No throat is regularized at c ≤ −1 ("Black universes with c ≤ −1, or equivalently 3πm ≥ 2b, have no throat in the
  static region"); V_eff is used as it stands.
- Time domain at c = −1 (Fig. 5, right, the curve for c = −1): "From Fig. 5 one can see that for c = −1 the effective
  potential has a negative gap, however, time-domain integration proves that in this case the black universe is stable."
  The method is cited as that of Gundlach et al. (§V.B, for the wormholes). **No initial data are described** in the
  text (full text read). The curve falls to about 0.03 by t ≈ 5 and stays within a factor of about 2 of that to t = 40
  (by eye, log axis; no numbers in the text). [Bro18]'s Fig. 5 shows time-domain profiles from [BKZ12] and says of the
  stable model: "the time-domain plot is asymptotically constant"; #232's 1D time domain gives a slope −0.037 over [45,
  60].
- Not in the paper: a count of growing modes, any growth rate, polar ℓ ≥ 1, the region beyond the horizon, the nonlinear
  problem ("One can mention a full nonlinear analysis of perturbations in all relevant cases.").

**Their horizon condition against #212's and S1's.**

| | [BKZ12] | #212 | S1 (#248, #232) |
|---|---|---|---|
| gauge | δR = 0, Eq. (35) | δφ = 0 | none: the invariant `Y = δr − (ρR/√2)δφ` |
| variable | Ψ; printed δφ = Ψe^{μ₃} (38), derived δφ = Ψ/R; `Ψ = −√2 u` (#212 S4, S8) | `u = Y/ρ` | Y |
| horizon condition | Eq. (43): ingoing, Ψ ∝ e^{Wr*}, finite at the horizon | ingoing: `Y ∼ ρ^{1+W/A₁}`, so `u → 0` | regular in Kruskal coordinates: at W = κ the branch `ρ^{1−W/A₁} = ρ^{1/2}`, `Y_h ≠ 0`, `u ∼ 1/ρ` |
| class admitted at c = −1 | Ψ finite at the horizon | `Y = O(ρ^{3/2})`, so `Y_h = 0`, c₁ = c₂ = 0 | c₂ free |
| the Y_h ≠ 0 solution seen in δR = 0 | `δφ = −φ′Y/R′ = −√2 Y/(Rρ)` (#212 §2 gauge map, with φ′ = √2/R², R′ = ρ/R), for `Y_ref ∼ ρ^{1/2}`: δφ ∼ ρ^{−1/2}, not finite | excluded | admitted |

So the horizon condition (43) on Ψ is #212's condition `u → 0` (the factor 1/R is 1 at the horizon, so the printed and
the derived (38) give the same condition). At c = −1, R′ vanishes at the horizon, so the gauge map from Y to δφ is
singular there. The last row of the table (the Y_h ≠ 0 solution seen in δR = 0) is by hand; the horizon exponents are
#232's D3a and the independent check of its review.

**What the same group says about the gauge.** These are statements for throats; that they apply at c = −1, where R′ = 0
on the horizon, is my reading and none of the papers says it.
- [BFZ11] §2.5 (the same first author as [BKZ12], eight months earlier): "The gauge δβ = 0 (the same as δr = 0) is
  suitable for describing the perturbations at any points except those where r′ = 0: these are throats and other
  critical points of r(u)." and "putting δr = 0, we forbid perturbations of the throat radius, while there is no
  physical reason for that". §2.6: "with nonzero δr on the throat. So it is necessary to pass on to another gauge", and
  Eq. (50): δr is finite in the gauge δφ = 0 when φ′ ≠ 0.
- [Bro18] §4.1: "The gauge δβ = 0 is convenient for calculations but looks doubtful when applied to configurations with
  throats." and "we suppose that the throat radius is not subject to perturbations, whereas its changes are in general
  admissible". The same review gives the boundary requirement (§4.3): "the perturbation δφ should remain finite in the
  whole space and also vanish at infinities or horizons". [Bro19] §3.5: "if we adopt for the horizons the same boundary
  condition δψ = 0 as is used at infinity (which looks reasonable)".
- [BIK26] §III.A (Nov 2026): "The gauge δβ = 0 is attainable only for perturbations with δβ = 0 at the throat, whereas
  the gauge δψ = 0 is always attainable. This is the source of all subtleties discussed below." In the gauge δψ = 0
  (their B, our Y), "B is regular everywhere".

None of them applies this to the horizon at the minimum of the area function, where [BKZ12] work in the gauge δR = 0.

**Restatements.** Eleven papers repeat the verdict without a new computation: [Bro18] (Fig. 5 caption), [BBS20] §5.3
("among the black-universe solutions there is a stable subset, in which the horizon coincides with the sphere of minimum
radius"; the boundary condition at a horizon there is "no reason to require there anything more than finiteness of δφ"),
1908.09126, 2009.12057, 2304.07383, 2305.00058, 2310.19037, 2510.18647, 2602.16972, 2604.20661 and 1308.5959 ("when the
horizon coincides with the throat, since it was shown in Ref. [21] that the solution is stable"). 2602.16972 and
2604.20661 treat massive test fields on the family, with the horizon at the throat for c = −1 (the critical scalar
charge); a test field is not the coupled radial perturbation. No paper of the 424 re-examines c = −1 with `Y ≠ 0` on the
horizon.

## 4. Horizon instabilities: the contrast for S3

| paper | system | horizon | what does not decay or grows | rate | data |
|---|---|---|---|---|---|
| [Are11a], [Are11b] | linear scalar, extreme RN | κ = 0 ("a black hole is called extreme if the surface gravity vanishes on H⁺") | transversal derivatives; higher derivatives blow up | power law in v | generic regular data on a slice crossing H⁺ |
| [Are15] | scalar wave equation on axisymmetric extremal horizons (assumption A2, Eq. 2: ∇_V V = 0) | κ = 0 | translation-invariant derivatives, from a conservation law on H⁺ | power law | "The genericity here refers to a condition on a section of H+" |
| [LMRT13] | massless, massive scalar and coupled linearized gravito-electromagnetic, extreme RN | κ = 0 | transverse derivatives | power law | "instabilities are present for initial perturbations supported outside the horizon, e.g. an ingoing wavepacket" |
| [MRT13] | nonlinear, spherical, Einstein–Maxwell–scalar | extreme RN | the scalar; generically the endpoint is a non-extreme RN | | |
| [CDLPV18] §3 | mass function of regular black holes, ingoing tail and outgoing flux | inner horizon, κ₋ | "the mass parameter in the region A grows exponentially, on a timescale determined by the surface gravity of the inner horizon" | e^{|κ₋|v} | generic, with ingoing tails |
| #248 | linear gravity + phantom scalar, c = −1 | non-extremal, κ = 2/π, at the area minimum | δr on H⁺ | `c₁ + c₂e^{κv}` | data that cross the horizon; c₂ = 0 for data supported away from the crossing (class R) |

The papers on κ = 0 horizons (also 1809.10037, 1606.08505, 1911.11164, 2112.09832, 2604.15260, 2610.11664; 2112.09832
uses the areal radius as the coordinate, so r′ = 1 on the horizon) assume extremality by definition; [MRT13] §1 gives
the reason ("for a non-extreme black hole, outgoing radiation at the horizon decays. However, for an extreme black hole,
κ = 0 so the horizon redshift effect is absent"). They concern a test scalar or the coupled perturbations of an extremal
horizon. The growth in S3 is of δr, at a non-extremal horizon, with the surface gravity as the rate; the data that
excite it differ from [LMRT13]'s (supported away from the horizon, there). No paper read has a non-extremal horizon with
an Aretakis-type or a κ-rate growth of the area radius on the outer horizon. The inner-horizon growth at κ₋ is a
different mechanism (blueshift of an ingoing tail against an outgoing flux, in the mass function).

## 5. Gauge fixing of the horizon: [HW13] and [Mars12] against the lemma

[HW13] §2.1 fixes the gauge for a perturbed stationary black hole in the vacuum, so that the surface at the old
horizon's coordinates is the event horizon to first order: "the linearized Raychaudhuri equation implies that, to first
order in λ, the expansion must be constant along the horizon. Thus, the perturbed expansion δϑ must vanish on the
horizon at all times". Then: "We now wish to argue that this condition can always be imposed." The argument displaces
the cross-section along the transverse null direction by f; the expansion changes by `C(f)` (Eq. 8); with the Killing
field `K = κ[u∂_u − r∂_r]` (Eq. 9) and the argument "if λ1 = 0, we would have ϑ(1) = 0, i.e. B1 would be a MOTS lying
outside the black hole, which is impossible" they obtain "λ1 is strictly positive", so that "we can uniquely solve the
equation C(f) = j for any smooth j". The statement is in Lemma 1 for "a stationary black hole spacetime satisfying the
vacuum Einstein equations". They add: "For another argument leading to the same conclusion in the more general context
of “non-evolving horizons”, see Prop. 3 of [26]", which is [Mars12].

[Mars12] (abstract): "We link these stability properties with the surface gravity of the horizon and/or to the existence
of minimal sections." The stability operator of a section S₀ of a totally geodesic null horizon is Eq. (11), `L ψ = −Δψ
+ 2 s_ℓ(∇ψ) + ψ(2 div s_ℓ − ½κ_ℓ θ_k + ½⟨k, tr N⟩)`; for a Killing horizon N = 0 ("This definition includes all Killing
horizons"). Theorem 1: "If H admits a minimal section (i.e. a section S0 with vanishing mean curvature vector H|S0 = 0),
then H is marginally stable", with κ_ℓ ≠ 0 constant; in the introduction: "non-evolving horizons with non-zero and
constant surface gravity are marginally stable if and only if they admit a section with vanishing total null expansion".

How these map on the lemma (by hand; not checked in code):

| lemma (#261) | literature |
|---|---|
| δr_VV = 0 on U = 0; δr = c₁ + c₂V; in the advanced time V = e^{κv}, so δr = c₁ + c₂e^{κv} | [HW13]: δϑ is constant along the horizon in the affine u (so δr = c₁ + c₂u); the flow of K maps u to ue^{κα} (Eq. 9 and the text before it), so u = e^{κv} up to a constant. The last step is the lemma's |
| c₂ shifts by −η(0)R′(ρ_h)ρ′(0) under U → U + η | [HW13] Eq. (8): `d/ds ϑ = C(f)`; [Mars12] Eq. (11): for a spherical section of a Killing horizon `L = −½κ_ℓθ_k`, so a deformation changes ϑ_ℓ by an amount ∝ θ_k ∝ R′(ρ_h) |
| R′(ρ_h) = 0: c₂ gauge invariant | [Mars12] Thm. 1: minimal section (θ_k = 0) ⇒ principal eigenvalue 0 (marginal stability); then `C(f) = j` cannot be solved for j = −δϑ|_B ≠ 0, the gauge step of [HW13] fails |
| the lemma holds for matter with T_VV = 0 on H⁺ | [HW13] vacuum, with the MOTS argument; [Mars12] Eq. (13) and Corollary 1 carry the Einstein tensor, with no energy condition for Prop. 2 |

The step from a gauge transformation of the metric perturbation to a deformation of the section is [HW13]'s; the
vanishing of L at a minimal section is [Mars12]'s; neither mentions δr, c₂, e^{κv}, a scalar field with a potential, or
a horizon at an extremum of the areal radius as a special case. [HW13]'s exclusion "MOTS lying outside the black hole,
which is impossible" is the standard statement that marginally trapped surfaces lie inside the black hole, which is
proved under the null energy condition; the c = −1 member violates that condition on either side of the horizon ([Bro18]
§2.2). That reading is mine, not [HW13]'s.

## 6. Other solutions with the minimum of the areal radius on a horizon

| solution | simple or extremal horizon | stability analysis in the paper | where |
|---|---|---|---|
| c = −1 black universe, ρ₀ ≠ 0 | simple, κ = 2/π | [BKZ12]: stable (the verdict above) | [BF06] after Eq. (14), [BKZ12] §II |
| ρ₀ = 0, c = −1 of the same family | double ("at c = −1 they merge into a double horizon") | none | [BF06] after Eq. (14) |
| the conformal image of the black universe, hybrid metric-Palatini | "a conformally mapped black universe with a horizon at the minimum radius" | by [BKZ12]'s result, no new computation | [BBS20] §5.3, §6 |
| the one-way wormholes of hybrid metric-Palatini gravity, charged, phantom sector | "either an EH or an extremal EH" at the throat; which, per solution, not read | "guiding future research into the stability properties" | [RJLR25] §V.C, §VI |
| examples 2 and 3 of generalized hybrid metric-Palatini gravity | "the horizon occurs at x = 0 and coincides with the minimum of r(x)" (m ≈ 0.212, 0.141); simple or double not stated in the caption | none | [BBS21] Figs. 4, 6 |
| Simpson–Visser, a = 2M | extremal ("a regular extremal black hole with a horizon at minimum radius") | [BIK26]: "the extremal case a = 2M is not analyzed separately"; "the regularity of perturbations at the degenerate horizon of the extremal configuration a = 2M remains to be established" | [SV19], [LRSSV21] §IV.A, [BIK26] abstract, §II.C, §VI |
| T-duality black bounce, `l₀ = 6M`; rotating Simpson–Visser at `l = r₊`; k-gravity black bounce `a = µ/b` | extremal for the first two ("degenerate (extremal horizon)"); not read for the third | none for the radial sector | 2603.05543 §II.C, 2201.01650, 2409.12101 §VI |

`horizon_lemma.md` already lists the c = −1 member and Simpson–Visser at a = 2M. The extremal ones are outside the lemma
by its item 4, which fits the Simpson–Visser metric: A = 1 − 2M/R is a function of the areal radius alone, so R′ = 0 on
a horizon forces κ = 0.

On the Simpson–Visser black-hole branch a < 2M, [BIK26] find a leading growing radial mode at every a and a growth rate
that tends, from both sides of a = 2M, to M·0.10121 (abstract: "the leading unstable black-hole and wormhole modes
approach the same limiting growth rate"); κ → 0 there. Here the limit of the rate is not the surface gravity, which is
the c = −1 family's behavior (#212). The paper says that this "does not establish the regularity of the perturbations at
the extremal horizon, which requires a separate analysis" (§V.D).

Hayward's view of a throat ([Hay09] §I): "a Morris-Thorne wormhole throat is a double trapping horizon, which will
generally bifurcate under a dynamic perturbation, such as someone crossing it", and §IV: "In static cases, a wormhole
mouth as defined above must be a double trapping horizon". The c = −1 horizon is a null hypersurface with r_U = r_V = 0
on it, a double marginal locus in the same sense. [Hay09] speaks of timelike throats of Morris–Thorne wormholes and of
generic bifurcation, not of a growth rate. [XCY25] (nonlinear double-null evolution of the Ellis wormhole) report that a
pulse "splits the throat into two distinct apparent horizons".

## 7. S4

No paper of the 424 gives a closed-form solution of a black-universe perturbation equation. [BKZ12] §V.A: "For nonzero
scalar potentials V (φ), analytical expressions for static perturbations are unknown. Therefore, such suitable solutions
must be found numerically under proper boundary conditions." Its Eq. (46) is the general static solution for V ≡ 0 only.
The screening pattern for closed forms (`exact_mode`) matched 7 papers: [BKZ12] (that sentence) and 6 on other systems
(test-field QNMs, Fisher solutions, a k-essence equation in 1908.09126). [BIK26] §VI recovers the Ellis rate `κa =
1.18197` (#212: 1.18197); no other rate that #212 gives is quoted in any of the 424. At c = 0 the S4 solution `Y = R` is
the variation of the throat radius b (#262 S2); [BFZ11] §2.6 and [BIK26] §III.A describe that mode in words (changing
throat radius).

## 8. Errors in the printed equations: what the literature does

- [BKZ12] printed (36)–(38) against #212's derivation (README §2 table): the coefficient of δφ′, δφ = Ψe^{μ₃}, two signs
  of (37). [Bro18], Eqs. (64)–(68), has: the δφ′ coefficient `−(2β′ + h′/h)` (for h = −1 phantom, `−2β′`, the sign of
  the derived form), `δφ = ψ e^{−η}, η′ = β′ + h′/2h` (for h = −1, 1/r, the derived form `Ψ/R`) and a U with `(V −
  e^{−2β})` and `+2φ′V_φ/β′` (the derived signs). Normalizations of V and φ differ from [BKZ12]'s; only signs and the
  factor e^{−β} are compared. By hand.
- Neither the INSPIRE record of [BKZ12] nor that of [BF06] has an erratum or related record (queries `q16`, `q17`); a
  title search for "erratum" with Bronnikov as author gives 0 hits.
- [BF06] Eq. (8): `κ = |A′(h)|/2 and A′(h) = (ρ0 − h)/r2 (h)`, while Eq. (7) `B′ ≡ (A/r²)′ = 2(ρ₀ − ρ)/r⁴` with A(h) = 0
  gives `A′(h) = 2(ρ₀ − h)/r²(h)`, as #212 S3 found. Checked again by hand from the printed (7) here. No later paper of
  the 424 repeats Eq. (8).
- The existing entries of `notes/literature.md` for [BF06] and [BKZ12] are consistent with the text re-read here; no
  change.

## 9. Sentences for Paper IV's `% NEEDS #260` comments

Wording follows the verdicts. Keys not in `main.tex` yet are given in `notes/literature.md`. `\ref{sec:prior}` stands
for wherever the coordinator puts the statement of what was searched.

1. Abstract, line 52 ("a sentence on what is new"): *Ref.~\cite{BKZ12} concludes that the $c=-1$ member is stable; that
   holds for the perturbations their gauge and horizon condition admit, which have $c_2=0$. We found no statement of the
   growing direction in the papers we searched (Sec.~\ref{sec:prior}).* No sentence says new; "found no statement" is
   the verdict "not found in the 424 papers". If a count belongs in the abstract: "in the 424 papers we searched".
2. Introduction, line 64 ("That class is this member."): *That class is this member ($c=-1$, $b=3\pi m/2$ in their
   Conclusion 2~\cite{BKZ12}; the horizon at the throat for $3\pi m=2b$ is already in~\cite{BF06}).*
3. Introduction, line 66, after "does not see.": *Their conclusion, stated for polar monopole perturbations (\S V C
   of~\cite{BKZ12}), is contradicted as worded. Their gauge $\delta R=0$ (Eq.~(35)) and their quasinormal condition on
   the horizon (Eq.~(43)) admit only $\Psi$ finite there (in our variable $u=Y/\rho$ finite, and $u\to0$ for an ingoing
   mode), which sets $Y=0$ on the horizon (Sec.~\ref{sec:cauchy}). The same gauge is called unsuitable at critical
   points of $r$ in~\cite{BFZ11}, for throats.* The last sentence is for throats in those papers; apply it to the
   horizon only as "the same reasoning", not as their statement.
4. Line 146 (the printed (36)–(38) disagree): *The printed Eqs.~(36)–(38) of~\cite{BKZ12} differ from the derivation in
   the sign of the coefficient of $\delta\varphi'$, in $\delta\varphi=\Psi e^{\mu_3}$ (we find $\delta\varphi=\Psi/R$)
   and in two signs of Eq.~(37); Eqs.~(40) and~(41) agree. Bronnikov's review~\cite{Bro18} prints Eqs.~(64)--(68) with
   the signs we derive. We found no erratum on INSPIRE. Eq.~(8) of~\cite{BF06} gives $A'(h)=(\rho_0-h)/r^2(h)$, half the
   value that its Eq.~(7) gives; no later paper that we read repeats it.*
5. Conclusions, "Prior work" item, line 591: append *Ref.~\cite{BKZ12} gives no initial data for its time-domain run at
   $c=-1$ (its Fig.~5), so we make no claim about the data of that run. The verdict is repeated without a new
   computation in~\cite{Bro18,BBS20} and in nine further papers.*
6. Conclusions, the mechanism item, line 594: replace the sentence on which other known solutions have such a horizon
   by: *In the papers we read, the minimum of the areal radius lies on a horizon in the $c=-1$ member~\cite{BF06,BKZ12},
   in its conformal image in hybrid metric-Palatini gravity~\cite{BBS20}, in the one-way wormholes of~\cite{RJLR25}, in
   two examples of~\cite{BBS21}, and in the extremal Simpson--Visser black bounce $a=2M$~\cite{SV19,LRSSV21}, which the
   lemma excludes; none of them is analyzed for the growth of Lemma~\ref{lem:horizon}, and~\cite{BIK26} lists the
   regularity of perturbations at the $a=2M$ horizon as open. In geometric language the lemma says that the horizon's
   sections are minimal, so that the stability operator of~\cite{Mars12} has zero principal eigenvalue (its Thm.~1) and
   the gauge condition $\delta\vartheta=0$ of~\cite{HW13} cannot be imposed; we found no paper that states this for
   $\delta r$ or that draws the growth at $\kappa$.* The geometric sentence credits S3's mechanism; the verdict is
   *contained in [Mars12] and [HW13]*, so the paper should not call the gauge invariance of $c_2$ new.

## What this does not show

- The patterns are fixed phrases; a paper that words S1–S4 differently from every pattern would be missed. The 59 read
  papers are the ones I chose from the citation lists and the searches, not every paper on the subject.
- The 52 citers of [BKZ12] without an arXiv id, the 2017 chapter on trapped ghosts, and the MOTS-stability papers named
  in Sec. 2 were not read.
- Whether [BKZ12]'s time-domain run contained data with c₂ ≠ 0 is unknown: the paper does not say. The reading that
  their class is `Y = 0` on the horizon rests on Eqs. (35), (43), on #212's identification Ψ = −√2u, and on #212's
  agreement with their verdict.
- The mapping in Sec. 5 is by hand: [Mars12] Eq. (11) was not rederived for the double-null chart, and the
  proportionality constants of "∝ R′(ρ_h)" were not compared.
- WebSearch returns at most ten links per query; its own summaries were not used.

## Recommended tasks

1. **The stability operator in the lemma's chart.** *What:* derive, for the spherical Killing horizon in #248's Kruskal
   chart with a scalar source, the first-order change of the outgoing expansion under U → U + η (that is, `L ψ` of
   [Mars12] Eq. (11), `−½κθ_k`), and compare with the lemma's `−η(0)R′(ρ_h)ρ′(0)`. *Why:* turns Sec. 5's by-hand mapping
   into *checked symbolically*. *Acceptance:* exact agreement at c = −1, −3/2 and Schwarzschild; a sign mutation in θ_k
   fails. *Depends on:* #261.
2. **Reproduce [BKZ12] Fig. 5's c = −1 curve and find its data class.** *What:* #232's characteristic 1D solver on Ψ
   with a Gaussian in Ψ on both rays, at c = −1; compare the shape of |Ψ| to t = 40 with the figure (digitized from the
   PDF). *Why:* the paper does not state its data; the curve shows the class if it matches. *Acceptance:* a data class
   that reproduces the dip at t ≈ 4 and the plateau to 20%, with the class's c₂ stated; a class with c₂ ≠ 0 shown to
   differ. *Depends on:* #232, #248.
3. **The extremal Simpson–Visser horizon.** *What:* the horizon sector at a = 2M in the same double-null scheme (κ = 0),
   against [BIK26]'s open item and the Aretakis-type growth. *Why:* the lemma excludes it, [BIK26] §VI lists it as open.
   *Acceptance:* the horizon law of δr at a = 2M (polynomial or constant in v), checked against an exact pure-gauge
   solution. *Depends on:* #223, #248.
4. **Read the MOTS-stability sources.** *What:* Andersson–Mars–Simon and Booth–Fairhurst, for any statement about
   marginal stability at minimal sections of horizons with NEC-violating matter. *Why:* closes the "not read" item for
   S3. *Acceptance:* a list of statements with section and equation, or "none". *Depends on:* #260.

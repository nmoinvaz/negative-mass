# Novelty check: are the candidate results already published?

Search run 2026-10-22. Raw log: `search_log.md` (same directory).

Tools used:
- WebSearch, about 50 queries in standard and extended modes.
- INSPIRE-HEP API citation trails, used as a stand-in for Google Scholar "cited by":
  - GHI hep-th/0409307: recid 660779, 61 citers
  - Gleiser–Dotti gr-qc/0604021: recid 713860, 59 citers
  - Dotti–Gleiser 0809.3615: recid 797174, 21 citers
  - Mbarek–Paranjape 1407.1457: recid 1305085, 26 citers
  - Shirokov 2608.24577: recid 3195733, 1 citer
  - Ivanov 1204.1269: recid 1104736, 3 citers
  - Novikov et al. 1807.06468: recid 1682851, 7 citers
- Every citer was screened by title. Plausible ones were opened via abstract, ar5iv/HTML, or PDF.

Note on INSPIRE: `refersto:arxiv:<id>` returns garbage for some IDs (about 25k unrelated hits). Only `refersto:recid:<n>` was trusted.

## Verdict table

| | Candidate | Verdict | Closest prior work |
|---|---|---|---|
| A1 | Exact zero mode φ₀ = 1/ψ₀ of the intertwined operator Ĥ | **Partially known** (implicit) | Dotti–Gleiser 2008 (0809.3615) §IV.1: the intertwiner maps 1/ψ to a solution of Ĥ; it is the SUSY-QM zero mode "missing in the spectrum" |
| A2 | Closed-form count #neg = [R<r_s] + [κ<κ_c(R)] for a truncated exterior with a Robin condition; thresholds r_s, r_* (l=2: 1.5, 2.279 \|M\|) | **Apparently new** | None found. Dotti 2022 review (2201.03640) §V *conjectures the opposite* |
| A3 | The ~l⁴ growth is a property of the point singularity; any finite star with l-independent κ caps it | **Apparently new** | Not stated anywhere found. DG08 has k ∝ l⁴ but no discussion of a cap |
| A4 | Non-radial stability of regular negative-mass stars | **Open in the literature** | Only radial results exist: Johnson–Paranjape et al. 2020 (radial, thin wall); de Aguiar Alves–Landulfo–Arderucio Costa 2025 (radial, barotropic ⇒ unstable) |
| B | Radial stability of phantom (ε=−1) boson stars; turning point holds in N, flips in M; first law dM = εω dN | **Apparently new** | No linear stability analysis of ε=−1 complex boson stars found. Shirokov 2608.24577 gives only nonlinear evolution of one star; Bernal et al. 2603.04400 (a different ghost model) disperses |
| C1 | Separation evolution d̈ = −E_tot/d²; rigidity \|μ₋\| = μ₊(1+μ₊/d); midpoint acceleration m/d² | **Known** (Ivanov 2012) | Ivanov, Phys. Usp. 55, 1232 (2012), arXiv:1204.1269: same equations by a different route |
| C2 | First-principles derivation from the signed-mass EIH Lagrangian; velocity-dependent 1PN formulas; Lorentz-covariance check; signed conserved integrals; application to Shirokov | **Apparently new**, as a derivation and an application | Ivanov explicitly leaves the EIH treatment "for possible future work" |
| D | Witten-type bounds with a signed source | **Not novel** (unchanged) | Shi–Tam 2004 (math/0406559); nothing newer in 2024–26 changes this |
| E | Responses to Shirokov 2608.24577 | **None exist** | Only a self-citation (Shirokov & Nachevsky 2610.13206). Earlier NR work on opposite-sign masses: Cao 2012, IJMPD 21, 1250061 |

---

## A. Negative-mass Schwarzschild: zero mode, truncated count, thresholds, l⁴ cap

**Queries used (brief).**
- "negative mass Schwarzschild perturbations stability boundary condition star regular interior" (extended)
- "Gleiser Dotti negative mass … later paper"
- "negative mass star stability non-radial perturbations"
- "Zerilli intertwining … zero mode 1/psi … Sturm"
- "negative mass thin shell/bubble stability Zerilli"
- "negative mass exterior matched interior star … unstable mode survives regular source" (extended)
- "S-deformation negative mass naked singularity"
- "Darboux transformation Schwarzschild 1/psi partner"
- INSPIRE trails of GHI, GD06, DG08 and Mbarek–Paranjape (about 165 citing records in total, screened by title)
- INSPIRE title search "negative mass" (94 records)

**Closest prior work.**

1. **Dotti & Gleiser 2008**, CQG 26, 215002 (2009), arXiv:0809.3615 (https://arxiv.org/abs/0809.3615).
   - Sec. IV.1, SUSY-QM remark: the partner solution "ψ^ωo(j=2) := 1/ψ_ωo(j=1)" satisfies Ĥ(1/ψ) = ω₀²/ψ. At ω₀ = 0 this is our φ₀. They note it is "missing in the spectrum" of Ĥ.
   - They never use it to count modes. The count of exactly one unstable mode is "numerical evidence".
   - They impose the boundary condition only at r = 0. Their artificial Robin condition at r = 10⁻⁴ is a numerical device.
   - **Conclusion: the φ₀ = 1/ψ₀ identity is implicit in DG08 and must be credited to them.** The node, critical-point and Sturm consequences (r_s, r_*, κ_c(R) and the count formula) do not appear.
2. **Dotti 2022**, "Linear Stability of Black Holes and Naked Singularities", Universe 8, 38, arXiv:2201.03640 (https://arxiv.org/abs/2201.03640).
   - §III.A covers M < 0. §V (Conclusions) states: "If the singularity were replaced by a small spherically symmetric negative mass distribution, away from r=0 … this potential would not change substantially, and an unstable mode could be found anyway."
   - It also says: "An analysis of these cases suggests that the singularity plays no role in the instability."
   - **This is the statement our result A2/A3 qualifies.** The instability survives only when R < r_*(l) (|M|/R > 0.439 for l = 2), or when the surface condition is softer than κ_c ≈ −(l+1)/R. For dilute stars and κ ≥ 0 it is absent, and the unbounded-in-l growth disappears. The paper gives no calculation behind the suggestion.
3. **Gleiser & Dotti 2006**, CQG 23, 5063, arXiv:gr-qc/0604021 (https://arxiv.org/abs/gr-qc/0604021).
   - Exact unstable modes; the Zerilli pole at r = −3M/λ (our r_s); static zero modes considered "not small".
   - No regular interior is discussed.
4. **Cardoso & Cavaglià 2006**, PRD 74, 024027, arXiv:gr-qc/0604101 (https://arxiv.org/abs/gr-qc/0604101).
   - Identifies the unstable modes as algebraically special, ω = i(l−1)l(l+1)(l+2)/(12|M|).
   - Negative-mass SdS is unstable; for SAdS the answer depends on boundary conditions.
   - No matter interior and no finite-radius cut.
5. **Dotti, Gleiser, Pullin, Ranea-Sandoval, Vucetich 2008 review**, arXiv:0810.0025: summarizes; nothing on stars.
6. **Gibbons, Hartnoll & Ishibashi 2004**, PTP 113, 963, hep-th/0409307: a one-parameter family of boundary conditions at the singularity. They do not tie the boundary condition to any resolution by matter.
7. **Giveon, Kol, Ori & Sever 2004**, JHEP 0408:014, hep-th/0401209.
   - Resolution of the negative-mass Schwarzschild singularity for *scalar* waves by smoothing. Negative-mass Schwarzschild is "fully reflecting".
   - Same philosophy as our observation that "the DG extension is the R → 0 limit of a regular cut-off". Cite as precedent for the scalar case.
8. **Regular negative-mass stars: radial only.**
   - Mbarek & Paranjape 2014, PRD 90, 101502, arXiv:1407.1457: existence; DEC must fail.
   - Johnson, Paranjape, Savard & Tapia-Arellano 2020, GRG 52, arXiv:1910.01774: thin-wall bubbles in dS, stable against *radial wall motion only*.
   - Novikov, Bisnovatyi-Kogan & Novikov 2018, PRD 98, 063528, arXiv:1807.06468: negative-density, positive-pressure stars. The paper says: "We did not consider here the issues of the stability."
   - **de Aguiar Alves, Landulfo & Arderucio Costa 2025**, PRD 111, 044027, arXiv:2408.00154: "any barotropic negative-mass star must be dynamically unstable". The analysis is radial (Chandrasekhar variational); non-radial instability is only heuristic. It does not cite Gleiser–Dotti.
   - Paranjape 2025, "Puzzles with Negative Mass" (Springer proceedings, ISBN 9783031912665, pp. 167–173): calls dynamical stability of the smoothed star "possible", i.e. open.
9. **Methods precedent** (cite for the technique, not the result): Araneda & Dotti 2017, PRD 96, 104020, arXiv:1611.03534, on Robin conditions at a timelike boundary producing a single unstable state; S-deformation literature (Kimura 2017, CQG, arXiv:1706.01447; Kodama–Ishibashi).

**Verdict.**
- A1: partially known; credit DG08.
- A2 and A3: apparently new, and they directly address the open claim in Dotti 2022 §V.
- A4: no non-radial stability analysis of any regular negative-mass star was found.

**Recommended framing.**
- Present φ₀ = 1/ψ₀ as "the zero-energy partner solution already implicit in Dotti–Gleiser's intertwiner [0809.3615, §IV]". Claim as new the *use* of it: the node at r_s, the critical point r_*, the closed form κ_c, and the exact count for a truncated Robin problem.
- State explicitly that this quantifies, and partly contradicts, the suggestion in Dotti 2022 §V that a regular negative-mass source would still be unstable.
- Keep the caveat that κ is a stand-in for the matter. The physical κ_eff(l, ω) of a real star is not computed.
- Must cite: GHI 2004, GD 2006, DG 2008, Cardoso–Cavaglià 2006, Dotti 2022, Giveon et al. 2004, Mbarek–Paranjape 2014, Johnson et al. 2020, de Aguiar Alves et al. 2025.

## B. Radial stability of phantom (ε = −1) boson stars

**Queries used (brief).**
- "phantom boson stars negative ADM mass stability" (extended)
- "Dzhunushaliev Folomeev phantom boson stars nontrivial topology PRD 90 124038"
- "Compact and extended objects from self-interacting phantom fields"
- "ghost complex scalar boson star negative energy … Q-ball" (extended)
- "boson star negative gravitational coupling turning point" (extended)
- "negative Newton constant boson stars Q-balls"
- "negative mass boson star" (extended)
- "gravitationally repulsive boson star negative active mass" (extended)
- "Kleihaus Kunz boson star negative ADM mass phantom"
- "Q-ball stability gravity … dN/dω"
- "catastrophe theory boson stars"
- INSPIRE title and abstract searches ("phantom" + "boson star", "negative mass" + "boson star" or "soliton"): 0 hits. INSPIRE's abstract index may be incomplete.

**Closest prior work.**
1. **Shirokov 2026**, arXiv:2608.24577 (https://arxiv.org/abs/2608.24577).
   - Constructs the ε = −1 star and claims it is "the first asymptotically flat body of negative ADM mass evolved in numerical relativity".
   - The lone star survives to t = 1000 with a peak decay of 0.63 % per 100 time units. He says this is not asymptotic stability.
   - No linear analysis, no N, no first law, no turning point. "Stability against sideways perturbations" is listed as open.
2. **Dzhunushaliev, Folomeev, Hoffmann, Kleihaus & Kunz 2014**, PRD 90, 124038, arXiv:1409.6978.
   - A *canonical* boson field plus a *real* phantom field supporting a wormhole core. The Ellis-throat instability persists.
   - Not a negative-active-mass boson star.
3. **Dzhunushaliev, Folomeev, Myrzakulov & Singleton 2008**, JHEP 0807:094, arXiv:0805.3211.
   - Static *real* phantom solitons with overall negative mass; reported unstable.
   - This is the closest "negative-mass phantom soliton". It has no U(1) charge or frequency, so there are no N or ω turning points.
4. **Dzhunushaliev et al. 2016**, PRD 94, 024004, arXiv:1606.07304: phantom balls from two real phantom scalars.
5. **Bernal, Jaramillo, Montiel-Hernández, Núñez & Sanchis-Gual 2026**, arXiv:2603.04400.
   - A complex ghost scalar confined in neutron stars. Their "free phantom boson star" (Fig. 8) disperses.
   - Their model differs: an overall sign flip with V = −μ²|χ|² + …, so the KG equation is changed.
   - No N, no first law, no perturbative analysis. Cite to distinguish the model: ours, Shirokov's, has an unchanged KG equation and only the Einstein-equation coupling flipped.
   - Their predecessor: Jaramillo, Jiménez-Vázquez & Núñez, PRD 107, 064061 (2023), arXiv:2302.07889 (real ghost in a boson star).
6. **Turning-point / catastrophe literature** (all with ordinary G > 0):
   - Gleiser & Watkins 1989; Lee & Pang 1989; Kusmartsev, Mielke & Schunck 1991; Jetzer; Kain 2021.
   - Tamaki & Sakai 2011, arXiv:1105.1498: criterion in Q; κ = m²G/λ is assumed positive, and negative κ is never considered.
   - Sorkin 1982 (ApJ 257, 847) for the general turning-point theorem.

**Verdict: apparently new.** No paper was found that does linear radial stability of complex boson stars with a sign-flipped gravitational coupling, or states dM = εω dN, or notes that the turning-point rule is ε-invariant only in N. The first-law and turning-point observations are close to a corollary of standard proofs (Sorkin; Lee–Pang) once ε is kept. A referee may call them "expected". The numerical stability result for the campaign and strong-field stars is the solid new content.

**Recommended framing.**
- Present as "the first linear (radial) stability analysis of a negative-ADM-mass boson star".
- The first law and the N-form turning-point criterion follow from the standard derivations with ε kept. Credit Sorkin, Lee–Pang and Kusmartsev–Mielke–Schunck, and claim the ε = −1 numerics and the M-form sign flip.
- Contrast with de Aguiar Alves et al. 2025: their theorem needs a barotropic fluid, and scalar-field stars evade it. This consistency check is worth one sentence.
- Contrast with Dzhunushaliev et al. 2008, whose real phantom solitons are unstable.
- Must cite: Shirokov 2026, Dzhunushaliev et al. 2008/2014/2016, Bernal et al. 2026, de Aguiar Alves et al. 2025, the standard turning-point references, Liebling–Palenzuela (LRR).

## C. 1PN two-body dynamics with one negative mass

**Queries used (brief).**
- "negative mass runaway post-Newtonian two-body" (extended)
- "EIH negative mass dipole relativistic correction uniformly accelerated" (extended)
- "Price 1993 Negative mass can be positively amusing"
- "Bonnor Swaminarayan … mass relation" (extended)
- "Bondi 1957 … masses must differ slightly"
- "negative mass Einstein-Infeld-Hoffmann Lagrangian runaway" (extended)
- "opposite signs masses post-Newtonian Havas Goldberg" (extended)
- "Fermi-Walker positive negative mass pair proper acceleration"
- Read Forward 1990 (PDF), Bondi 1957 (PDF), Ivanov 2012 (full text), Bičák–Kofron 2009, the Pravda–Pravdová review, and Trivedi–Loeb 2026.

**Closest prior work.**
1. **Ivanov 2012**, "On relativistic motion of a pair of particles having opposite signs of masses", Phys. Usp. 55, 1232–1238 (2012) [UFN 182, 1319], arXiv:1204.1269 (https://arxiv.org/abs/1204.1269).
   - With m_{1,2} = ±(μ ± Δμ):
     - comoving energy E_c = 2Δμ + μ²/D (eq. 11);
     - D̈ = −2Δμ/D² − μ²/D³ in Fermi–Walker coordinates (eqs. 16, 21);
     - midpoint acceleration a = μ/D_FW² (eq. 19);
     - constant separation requires Δμ < 0 with |Δμ| = μ²/(2D). The text has a typo, "D_FW = 2|Δμ|".
   - **These are exactly our d̈ = −E_tot/d², E_tot = μ₊ + μ₋ − μ₊μ₋/d, the rigidity condition |μ₋| = μ₊(1 + μ₊/d), and a_mid = m/d².**
   - His derivation comes from conservation of the four-momentum norm plus Fermi–Walker kinematics. He calls the μ²/D³ term "technically … post-Newtonian" and says an EIH / Havas–Goldberg treatment is "left for a possible future work".
   - Only 3 citers: Novikov et al. 2018, arXiv:2411.00933, arXiv:1406.0841. **Shirokov does not cite it.**
2. **Forward 1990**, J. Propulsion Power 6, 28 (https://lweb.cfa.harvard.edu/~loeb/forward1990.pdf): citing Bonnor–Swaminarayan, "the negative mass should be very slightly larger than the positive mass, depending on the initial separation". Qualitative statement of the rigidity condition.
3. **Bondi 1957**, RMP 29, 423 (https://www.januscosmologicalmodel.com/pdf/bondi1957.pdf): "the masses are not quite equal and opposite … since their accelerations must be unequal in the uniformly accelerated model." Qualitative.
4. **Bičák & Kofroň 2009**, GRG 41, 153, arXiv:0901.3971: Newtonian limit of Bonnor–Swaminarayan and dipole solutions. Masses are equal and opposite at Newtonian order, and the paper explicitly does no post-Newtonian work.
5. **Bonnor & Swaminarayan 1964**, Z. Phys. 177, 240; Bičák, Hoenselaers & Schmidt 1983: exact solutions (already in `notes/literature.md`).
6. **Trivedi & Loeb 2026**, arXiv:2605.10976: Newtonian plus quadrupole only, no 1PN.
7. **Price 1993**, Am. J. Phys. 61, 216: paywalled, not read. It is a pedagogical Newtonian treatment and very unlikely to contain 1PN.
8. Petit & d'Agostini 2014, arXiv:1408.2451 (bimetric, Newtonian): not relevant except as a critique-camp reference.

**Verdict.**
- C1 is **known** (Ivanov 2012, with qualitative antecedents in Bondi 1957 and Forward 1990).
- C2 is **apparently new** as a first-principles result:
  - the signed-mass EIH derivation, with the argument that m → μ continuation is exact through 1PN including regularisation and effacement;
  - "no 1PN correction to the proper midpoint acceleration", checked to O(V²);
  - the velocity-dependent formulas;
  - the signed E, P, J and G integrals and the centre-of-energy argument for E_tot = 0;
  - the application to Shirokov: E_tot > 0 implies contraction, and 1PN cannot produce a·d²/M̄ ≈ 1.08.
- Ivanov's derivation is heuristic about the μ²/D³ term, so an EIH confirmation has some value. It is a modest addition.

**Recommended framing.**
- Do **not** claim the rigidity condition or the separation law as new. Write "we recover, from the EIH Lagrangian with signed masses, the relations obtained by Ivanov [Phys. Usp. 55, 1232 (2012)] from energy–momentum conservation", and credit Forward 1990 and Bondi 1957 for the qualitative statement.
- Claim the first-principles derivation and its use to diagnose Shirokov's runs.
- Shirokov's "mass ladder" flips at ratio 0.995 against a predicted 1.000. The 1PN / Ivanov prediction is 1 + μ₊/d ≈ 1.0014. Both are well inside his ±0.03 uncertainty, which is worth a footnote.

## D. Witten-type bounds (brief re-check, 2024–2026)

**Queries used.**
- "lower bound ADM mass negative energy density Witten spinor 2025 OR 2026" (extended)
- "ADM mass lower bound scalar curvature negative somewhere integral estimate 2025 2026 Shi Tam generalization" (extended)
- "Witten positive energy ghost/phantom two scalar fields opposite sign"
- "Bondi mass ghost scalar null infinity"
- "positive mass theorem negative energy density L^{3/2} Sobolev spinor 2026" (extended)

**Found, none of which changes the verdict:**
- Cecchini–Lesourd–Zeidler, IMRN (2024), arXiv:2307.05277 (already cited in phase 4).
- arXiv:2509.05244, positive energy for spin initial data with charge (assumes DEC).
- arXiv:2604.24746, hyperboloidal / spacetime PMT in all dimensions.
- arXiv:2508.17585, creased initial data.
- arXiv:2606.23529, low-codimension singularities.
- Caridi, Corelli & Pani 2026, arXiv:2603.13210, NR collapse of phantom scalar wave packets: always disperse, no trapped surfaces, no negative-mass endpoint. Relevant context for gap 6.

**Verdict: unchanged.** Lemma A is folklore and Lemma B is a spinor analogue of Shi–Tam 2004 (math/0406559). Do not claim novelty. The *refutation* (E₋ is not conserved and the runaway pair makes W₋ grow) and the 0.11 % tightness on the phantom star can be presented as worked examples. Add Caridi–Corelli–Pani 2026 to the citations.

## E. Shirokov 2608.24577: responses and other NR runaway simulations

**Queries used.**
- INSPIRE `refersto:recid:3195733`
- WebSearch "Shirokov Bondi dipole numerical relativity comment OR response OR critique" (extended)
- "self-accelerating positive negative mass pair numerical relativity simulation 2026" (extended)
- "\"2608.24577\""
- "Bondi dipole phantom boson star … comment 2026"
- "numerical relativity negative mass simulation … 2024 2025 2026"

**Found.**
- Only one citing record: **Shirokov & Nachevsky 2026**, arXiv:2610.13206, "Binary Ellis–Bronnikov wormholes in 3+1 numerical relativity" (a self-citation, ghost-scalar wormholes).
- No comments, replies or errata. arXiv shows only v1 (2026-08-25).
- **Earlier NR of opposite-sign masses:** Cao Zhou-Jian 2012, "Numerical study of the interaction between positive and negative mass objects in general relativity", IJMPD 21, 1250061, DOI 10.1142/S0218271812500617.
  - Uses BSSN with moving punctures and excises the negative-mass naked singularity, with a custom outer boundary condition.
  - Recovers the Newtonian "chase" in the weak-field limit.
  - **Shirokov does not cite it.** His "first … evolved in numerical relativity" claim concerns an asymptotically flat, *regular* negative-ADM-mass body, which is defensible. But a runaway-pair NR study predates him: Cao's is singular and excised, and has 4 citations.
- Other related NR work: Caridi–Corelli–Pani 2026 (phantom collapse); the wormhole NR line cited by Shirokov (arXiv:2604.00071).

**Verdict.** Our code review would be the first public response to Shirokov.

**Recommended framing.** Cite Cao 2012 as prior NR of a ± mass pair, singular and excised. Cite Ivanov 2012 as an analytic prediction Shirokov missed, which agrees with our 1PN numbers.

---

## Limits of this search

- **arXiv full-text and listing API** (export.arxiv.org) was blocked by a permission timeout, so systematic arXiv abstract queries could not be run. The direct curl route to arXiv and INSPIRE was refused by the proxy; WebFetch to INSPIRE worked.
- INSPIRE abstract-field queries returned 0 hits even for terms that should match, so that index is unreliable. The citation trails (recid-based) are reliable.
- No Google Scholar access, so non-INSPIRE literature (pure mathematics, AJP, J. Propulsion) is covered only by web search.
- **Paywalled or not read:**
  - Price 1993 (AJP 61, 216);
  - Bonnor & Swaminarayan 1964;
  - Bičák, Hoenselaers & Schmidt 1983;
  - Cao 2012 full text (abstract only);
  - Paranjape 2025 proceedings (abstract only);
  - the full text of Dotti–Gleiser 2008 §V (ar5iv was truncated, so the conclusion section was not read).
- de Aguiar Alves et al. 2025 and Ivanov 2012 were read through tool summaries of the PDFs. The equations quoted for Ivanov were checked for internal consistency (one typo found in his text).
- Russian-language follow-ups to Ivanov (UFN) and Chinese-language follow-ups to Cao may exist outside INSPIRE.
- Results posted after about 2026-10-21 are not covered.
- Some journal references were filled in from memory and not checked against the publisher: DG08 CQG 26, 215002; GD06 CQG 23, 5063; Cardoso–Cavaglià PRD 74, 024027; Dzhunushaliev et al. JHEP 0807:094; Bičák–Kofroň GRG 41, 153; Mbarek–Paranjape PRD 90, 101502. Check them before submission. The arXiv IDs were all verified.

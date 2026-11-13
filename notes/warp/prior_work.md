# Prior work on the warp-drive results of #57 and #58

Issue [#60](https://github.com/nmoinvaz/negative-mass/issues/60), searched and read on 2026-11-13. Each work below is
recorded, with what it does and does not contain, in [`../literature.md`](../literature.md) under "Warp drives". This
note compares four results with that record: no scalar source and no perfect-fluid source
([`phantom_source.md`](phantom_source.md)), E₋ ∝ v and the flat-slice existence lemma
([`mass_bounds.md`](mass_bounds.md)).

## Verdicts

| result | where | verdict | closest prior work |
|---|---|---|---|
| No scalar field of either sign sources Alcubierre where f′ ≠ 0 | phantom_source.md, Answer and Lemma | partly contained, by implication | [Rod26]: Type IV on the whole axis; [Le26]: Type IV over most of the sampled wall; both numerical |
| No perfect fluid, timelike or null flow, sources it there | phantom_source.md, after the Lemma | partly contained | [SAR21a], [SAR21b], [ABV23]: u = n only; [Rod26], [Le26]: Type IV points |
| Minimal sector split needs E₋ linear in v | mass_bounds.md items 3–4 | partly contained | [LV04a], [LV04b]: NEC violation at O(v); [Le26]: exact v-scalings, deficit linear in \|j\| |
| Witten spinor exists on a flat slice when κ_K < 1 | mass_bounds.md item 5 | partly contained | [FNS11]: only tr K enters; [FL22]: the same Sobolev–Hölder argument for Dirac zero modes, with a sharper constant |

No work read contains any of the four in full, so the stop condition of #60 did not fire.

## 1. No scalar source

The result. For every C² profile and both signs ε, the Ricci tensor of the Alcubierre metric is not of the form
8π(ε ∂φ∂φ + V g) at any point where f′(r_s) ≠ 0, because of the identities (I1), (I2) and a rank lemma.

What prior work has.
- No paper was found that asks whether a scalar field, normal or phantom, can source the Alcubierre metric. The
  scalar-field papers on the trail are about wormholes, quantum fields on the warp background, or test fluids
  ([GD07], phantom accretion in two dimensions).
- A minimally coupled scalar has a stress tensor of Hawking–Ellis Type I (timelike or spacelike gradient) or Type II
  (null gradient). [Rod26] finds the Alcubierre G^a_b of Type IV along the whole axis of motion and over the equatorial
  cut, except an annulus of width ≈ 2Δ about the wall, for one parameter set. [Le26] finds the sampled Alcubierre wall
  Type I at 1.3 % of points at v = 0.5 and at 24.9 % at v = 2.5, the rest Type IV, and gives the criterion
  2|j| > |ρ_n + S_∥|. At a Type IV point no scalar can be the source, so these maps exclude scalars at those points.
  Neither paper draws that conclusion.
- On the axis ρ_n = 0 and S_∥ = 3ρ_n = 0, so [Le26]'s criterion holds wherever j ≠ 0. That is the complex pair
  −R_00 ± iR_01 of phantom_source.md, which [Rod26] reports numerically ("Type IV for all radii" on the axis).

Verdict: **partly contained, by implication**. The Type IV maps exclude a scalar source at most wall points, for sampled
parameters and numerically. The exclusion at every wall point, Type I points included, for every profile and both
signs, was not found.

## 2. No perfect-fluid source

The result. The rank lemma excludes R_ab = λ g_ab + s u_a u_b at every wall point, for any equation of state and any
timelike or null u. phantom_source.md calls this a corollary not yet compared with the literature.

What prior work has.
- [SAR20] (dust), [SAR21a] and [SAR21b] (perfect fluid, without and with Λ) and [ABV23] (dust, isotropic fluid,
  anisotropic fluid with two equal pressures) all fix **u = n**, the Eulerian normal, and allow a general β(t, x, y, z).
  The field equations then leave only vacuum (β independent of y and z, inviscid Burgers) or a non-localized β that is
  linear in y with μ ≤ 0, or complex. So no warp bubble has a perfect-fluid or dust source at rest in the Eulerian
  frame. Their key step is p = 3μ ([SAR21a] Eq. 3.4) or p_x = 3ρ ([ABV23] Eqs. 30–31), which is (I2) seen by a fluid
  with u = n (see "The identities" below).
- [BD20] fits two- and three-fluid mixtures and says one fluid gives an overdetermined system, without a proof.
- [CDK24] use the full anisotropic G_ab/8π, and leave a fluid with a general rest frame for later work.
- The Type IV maps of [Rod26] and [Le26] exclude every perfect fluid, with any u, at the Type IV points (a perfect fluid
  is Type I, or Type II for a null flow).

Verdict: **partly contained**. Contained for u = n and every profile ([SAR21a], [SAR21b], [ABV23]); contained for any u
at the Type IV points, numerically, for sampled parameters ([Rod26], [Le26]). A tilted fluid at a Type I wall point is
not covered by any work read. These works must be cited where the corollary is stated.

## 3. E₋ linear in v

The result. Over splits T = T₊ − T₋ with μ± ≥ |J±|, the pointwise minimum is μ₋ = max(0, −μ, (|J| − μ)/2), and its
integral E₋ grows like v at small v, driven by the wall's momentum density, while the Eulerian ∫μ dV grows like v².

What prior work has.
- The momentum density: [SSV22] Eqs. 4.7–4.8, f_i = (1/16π)(∇×∇×v)_i with v = −β, which is J_i of mass_bounds.md
  item 1 exactly. The Eulerian energy integral: [PF97] Eq. 26 and [LV04a], −(v²/12)∫r²f′²dr.
- Local energy-condition violation at O(v): [LV04a] Eqs. 15–19 (exact NEC along ± the direction of motion and
  transverse to it, with a term linear in v that dominates at small v), [LV04b] Eq. 8 (boosted observers see ρ < 0 at
  O(v)), repeated in [AL17]. [Le26] proves ρ_n, S_ij ∝ v² and |j| ∝ v exactly, and that the null witness deficit
  ρ_n + S_∥ − 2|j| is linear in |j|.
- Qualitative statements that a positive sector must outweigh the momentum flux: [Fu24] Sec. 5.2.
- DEC violation without a magnitude: [SSV22], and [BBV26] Corollary IV.37 by the positive energy theorem.

What it does not have. No work read splits T into two sectors that each satisfy μ ≥ |J|, minimizes the negative
sector, or integrates |J|. [LV04a] integrate the signed O(v) term, find zero, and quantify the net violation as
O(v²) ≈ −v²R²σ. That is a signed volume integral, not the integral of a pointwise minimum, so it does not conflict with
E₋ ∝ v. The two must not be conflated when [LV04a] is cited.

Verdict: **partly contained**. The O(v) local violation from the momentum density is prior work ([LV04a], [LV04b],
[SSV22], [Le26]). The minimal split, the integral E₋ and its comparison with the Witten weight were not found.

## 4. Flat-slice existence lemma

The result. On a flat slice ∇̃_i = ∂_i + ½K_ijγ^jγ⁰ and the Witten operator is D̃ = D + ½(tr K)γ⁰. With the sharp
Sobolev constant S = 3(π/2)^{4/3}, the map ψ ↦ −D⁻¹(½ tr K γ⁰(ε₀ + ψ)) is a contraction on Ḣ¹ when
κ_K = ‖tr K‖_{L³}/(2√S) < 1, which gives a unique Witten spinor and no decaying kernel. For tr K = 0, ε = ε₀.

What prior work has.
- That only tr K enters the Witten operator, in two-spinor form: [FNS11] Eq. 2.1. That constant spinors solve the
  Witten equation on maximal, intrinsically flat slices: [FNS11] Secs. 1 and 2.1.
- The analytic core: [FL22] recall as standard that a zero mode of γ·(−i∇ − A) needs ‖A‖²_{L³} ≥ S₃, by the Sobolev,
  Hölder and diamagnetic inequalities, which is the estimate behind κ_K < 1 with S₃ = S. They prove the sharp
  ‖A‖²_{L³} ≥ 3S₃ (Theorem 1), and Remark 4 extends it to any ψ ∈ L³ with |γ·∇ψ| ≤ |A||ψ|.
- E_ADM = 0 on flat slices and p_ADM = 0 for Alcubierre: [BBV26] Theorem IV.19.
- Existence of Witten spinors under the DEC (Parker–Taubes, Reula, Bartnik–Chruściel) is standard and was not re-read.
  No work read proves existence without the DEC under a smallness condition on K, or applies it to warp initial data.

Verdict: **partly contained**. The tr K structure and the maximal case are in [FNS11]. The smallness argument is the
standard one that [FL22] recall for Dirac zero modes. The statement for flat warp data, κ_K and v_c, was not found.

A possible sharpening, not checked here. A decaying kernel element ε = (a, b) of D̃ in mass_bounds.md's
representation gives σ·∇a = (i/2) tr K b and σ·∇b = −(i/2) tr K a. So ψ = a + ib solves σ·∇ψ = ½ (tr K) ψ and
χ = a − ib solves σ·∇χ = −½ (tr K) χ, and ψ and χ vanish together only if ε does. If they lie in L³, [FL22] Remark 4
with |A| = ½|tr K| excludes the kernel
for ‖tr K‖_{L³} < 2√(3S), a factor √3 above the contraction threshold. Multiplication by tr K ∈ L³ is compact from Ḣ¹
to L², so D̃ is Fredholm of index 0 and the absence of a kernel would give existence. The reduction, the L³ decay and
the Fredholm step are this note's reading, not checked. If they hold, each v_c of item 5 grows by √3 for existence,
while the contraction bound keeps κ_K < 1.

## Other items checked

- **The identities.** In the orthonormal frame of phantom_source.md, G_11 − 3G_00 = −(R_00 + R_11 + 2R_ρρ + 2R_ϕϕ)
  (*proved*, by hand, from G = R − ½Rη). With R_ϕϕ = 0 from (I1), (I2) is equivalent to G_11 = 3G_00. That is the
  identity G_zz = 3G_nn of [SSV22] Sec. 7.1, implicit in [LV04a] Eq. 16, and behind p = 3μ ([SAR21a]) and p_x = 3ρ
  ([ABV23]). R_ϕϕ = 0 itself was not found stated. [BD20] list no ϕ cross components, which is the rest of (I1).
- **E = P = 0** (mass_bounds.md item 2, *numerical evidence*) is [BBV26] Theorem IV.19 (iii). [LV04a], [SSV23],
  [Fu24] and [CDK24] state that the Alcubierre ADM mass vanishes, without computing it.
- **The Eulerian energy density** ρ_E (phantom_source.md) is Alcubierre's Eq. (19), and its integral is [PF97] Eq. 26.
- **The Witten identity at E = P = 0** (mass_bounds.md item 6): no prior computation found.

## Discrepancies

1. #60 names the dust and perfect-fluid authors "Santos-Pereira, Abreu and Rebouças". The third author is Marcelo B.
   Ribeiro (INSPIRE and the PDFs).
2. [FL22] would raise the existence threshold of mass_bounds.md item 5 by √3, if the reduction above holds. The note
   says "Beyond v_c this says nothing either way", which stays true of the contraction.
3. [LV04a] call the net NEC violation O(v²), while mass_bounds.md item 4 calls E₋ O(v). The quantities differ (§3);
   the notes should say so when they cite [LV04a].
4. [Rod26] find the equatorial wall annulus Type I, while [Le26] find the Alcubierre wall 98.7 % non-Type-I at
   v = 0.5. Their parameters and wall masks differ, and neither gives a closed form. An analytic type map from (I1)
   and (I2) would settle it (recommended task C in the #60 summary).
5. phantom_source.md says the perfect-fluid corollary "has not been compared with the literature", and mass_bounds.md
   says the prior work "has not been read". Both are now out of date. This task does not edit those files.

## Search log

All searches on 2026-11-13, from the coordinator's Mac.

INSPIRE-HEP API (`inspirehep.net/api/literature`):
- `refersto:recid:389445` (Alcubierre 1994): 368 citing records, all screened by title.
- `t warp`: 240 records, 102 with warp-drive titles screened. `a Santos-Pereira`: 130 records, the warp ones taken.
  `t "warp drive" and t scalar`: 0. `t "superluminal warp drive"`: 3.
- Full text: `"Alcubierre metric"` with `"phantom scalar"`: 0. `"warp drive"` with `"ghost scalar"`: 15, none a warp
  source. `"Alcubierre"` with `"Segre type"`: 7; Mattingly et al. 2020 (arXiv:2010.13693) leave the Segre type for
  later work. `"warp drive"` with `"rank one"`: 2 ([Le26] and arXiv:2606.22531, not about sources). `"warp drive"`
  with `"Witten spinor"`: 1, unrelated.
- Forward trails, `refersto:recid:` of [SAR20] 1811921, [SAR21a] 1843032, [SAR21b] 1912085, [ABV23] 2637024, [SSV22]
  1862524, [SSV23] 2089890, [LV04a] 652763, [BBV26] 3121358, [Rod26] 3094802, [Le26] 3122171, [Fu24] 2782713, Warp
  Factory 2774121 and [CDK24] 2794423: 325 records, 63 not already among Alcubierre's citers. None treats warp-drive
  sources, energy-condition splits or ADM quantities.
- Backward trails: the reference lists of [SSV22], [Rod26], [Le26] and [BBV26] (its "Errors" survey).

Web search (standard or extended):
- warp drive metric sourced by a scalar field; "warp drive" phantom scalar source; Alcubierre Hawking–Ellis type IV
  perfect fluid any four-velocity; minimal negative energy split into two DEC tensors; Alcubierre Ricci Segre type;
  Witten spinor existence without the DEC with small tr K; Frank–Loss zero modes and the L³ norm.
- Found [FL22] and [FNS11] this way. Nothing found on scalar sources of the Alcubierre metric, two-sector DEC splits,
  or Witten-spinor existence on warp data.

Reading. About 55 arXiv PDFs were downloaded and converted to text, and the sections named in the literature entries
were read. Google Scholar was not used. Results posted after 2026-11-13 are not covered.

## Not read, or read only in part

- [GD07] and González-Díaz, "Superluminal warp drive", Phys. Lett. B 653, 129 (2007): paywalled, no arXiv version,
  abstracts from INSPIRE only.
- [FL22]'s journal version: the arXiv version was read. Its journal reference comes from a web-search snippet.
- Read in abstract and setup only: [SAR] charged dust (arXiv:2102.05119), [ABV] spherical warps (arXiv:2305.03736),
  Lentz 2021 and 2023, Fell & Heisenberg 2021, Bobrick & Martire 2021, Buchert & Frackowiak 2026, Bolívar et al. 2026,
  Celmaster & Rubin 2025, and Warp Factory (grep only).
- Not read: Parker–Taubes 1982, Reula 1982 and Bartnik–Chruściel 2005 (existence under the DEC, standard); the
  Hawking–Ellis classification papers of Martín-Moruno and Visser; Olum 1998, Low 1999 and Coule 1998 (superluminal
  travel and energy); Hiscock 1997 (quantum effects); McMonigal et al. 2012 and Carneiro et al. 2022, whose PDF text did
  not extract.

## References

Full entries, with journals checked against INSPIRE, are in [`../literature.md`](../literature.md).

- [ABV23] Abellán, Bolívar, Vasilev, arXiv:2302.13826, GRG 55, 60 (2023).
- [AL17] Alcubierre, Lobo, arXiv:2103.05610, Fundam. Theor. Phys. 189, 257 (2017).
- [BBV26] Barzegar, Buchert, Vigneron, arXiv:2602.16495 (2026, preprint).
- [BD20] Béatrix-Drouhet, arXiv:2012.09941 (2020, preprint).
- [CDK24] Clough, Dietrich, Khan, arXiv:2406.02466, Open J. Astrophys. 7 (2024).
- [FL22] Frank, Loss, arXiv:2201.03610 (2022); JEMS 28, 305 (2026) *(unverified)*.
- [FNS11] Frauendiener, Nester, Szabados, arXiv:1105.5008, CQG 28, 185004 (2011).
- [Fu24] Fuchs et al., arXiv:2405.02709, CQG 41, 095013 (2024).
- [GD07] González-Díaz, Phys. Lett. B 657, 15 (2007), doi:10.1016/j.physletb.2007.09.041.
- [Le26] Le, arXiv:2602.18023, CQG 43, 215013 (2026).
- [LV04a] Lobo, Visser, arXiv:gr-qc/0406083, CQG 21, 5871 (2004).
- [LV04b] Lobo, Visser, arXiv:gr-qc/0412065 (2004, proceedings).
- [PF97] Pfenning, Ford, arXiv:gr-qc/9702026, CQG 14, 1743 (1997).
- [Rod26] Rodal, arXiv:2512.18008, GRG 58, 1 (2026).
- [SAR20] Santos-Pereira, Abreu, Ribeiro, arXiv:2008.06560, EPJC 80, 786 (2020).
- [SAR21a] Santos-Pereira, Abreu, Ribeiro, arXiv:2101.11467, EPJC 81, 133 (2021).
- [SAR21b] Santos-Pereira, Abreu, Ribeiro, arXiv:2108.10960, EPJ Plus 136, 902 (2021).
- [SSV22] Santiago, Schuster, Visser, arXiv:2105.03079, PRD 105, 064038 (2022).
- [SSV23] Schuster, Santiago, Visser, arXiv:2205.15950, GRG 55, 14 (2023).

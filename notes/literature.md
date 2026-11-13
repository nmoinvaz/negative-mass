# Literature notes

## Foundations
- **Bondi 1957**, Rev. Mod. Phys. 29, 423 — inertial / passive / active mass; runaway pair from Weyl metrics
  (negative-mass source cancels the conical defect). *Paywalled; not yet read in full.* Chapel Hill version:
  https://edition-open-sources.org/sources/5/24
- **Bonnor & Swaminarayan 1964**, Z. Phys. 177, 240 — exact uniformly accelerated pairs.
- **Bičák, Hoenselaers & Schmidt 1983**, Proc. R. Soc. A 390, 397 & 411 — solutions without nodal singularities
  (part II: self-accelerating particles). Only pairs with one negative mass move freely without struts.
- **Bonnor 1989**, GRG 21, 1143 — assumes no runaways ⇒ all mass negative; negative-mass cosmologies;
  footnote: Bondi's same-sign-inertial / opposite-active case violates the third law but is allowed in GR.
  https://januscosmologicalmodel.com/pdf/bonnor1989.pdf
- **Forward 1990**, J. Propulsion Power 6, 28 — runaway conserves momentum and energy.

## Linear stability of negative-mass Schwarzschild
- **Gibbons, Hartnoll & Ishibashi 2004**, hep-th/0409307 — vector modes stable; scalar modes: family of BCs
  q ∈ ℝℙ¹, q_C = 2 + log(98/9), Friedrichs extension "stable". Flaw: truncation at r_C (see phase0/RESULTS.md).
- **Gleiser & Dotti 2006**, gr-qc/0604021 — kinematic singularity at r = 3μ/(2λ); exact smooth unstable modes.
- **Dotti & Gleiser 2008**, arXiv:0809.3615 — intertwiner Ψ̂ = (∂ₓ − ψ₀'/ψ₀)Ψ_z; unique self-adjoint
  extension; instability excitable by generic data. Corollary: linearly unstable.
- Related: Dotti–Gleiser–Pullin gr-qc/0607052 (super-extremal RN/Kerr); Cardoso–Cavaglia gr-qc/0604101
  (negative-mass Schwarzschild–(A)dS).

## Regular negative-mass sources
- **Mbarek & Paranjape 2014**, arXiv:1407.1457 — any nonsingular matter approaching negative-mass Schwarzschild
  must violate the DEC; DEC-satisfying perfect-fluid bubbles exist in de Sitter. Resolves the apparent tension
  with Shirokov (his phantom violates the NEC).
- Belletête & Paranjape 2013, arXiv:1304.1566; Mbarek, Savard, Paranjape, arXiv:1910.01774 (stable thin-wall bubbles).

## Quantum
- **Cline, Jeon & Moore 2004**, hep-ph/0311312 — ghost vacuum decay via gravitons; divergent phase space;
  needs a Lorentz-violating cutoff. Carroll, Hoffman & Trodden 2003 (PRD 68, 023509).

## Numerical relativity
- **Shirokov, Aug 2026**, arXiv:2608.24577 (single-author preprint, not peer reviewed) — two complex scalars,
  phantom sector enters Einstein eqs with a minus sign (positive inertial/passive, negative active mass).
  Runaway survives full GR: a ∝ d^−2.03, signed momentum conserved, no detectable radiation; lone phantom
  star to t = 1000. Open: off-axis kicks, rigidity past t ≈ 200, relativistic-speed radiation; N = 192 test
  predicts d(400) ≈ 10.7.
- **Ghost-wormhole thresholds vs. Choptuik critical collapse** (Shinkai–Hayward gr-qc/0205041, González–Guzmán–Sarbach
  0806.0608/0806.1370, normal + ghost σ-models, phantom collapse): see `notes/critical_ghost.md`.

## Warp drives: matter sources, energy conditions, ADM quantities
Read for #60 on 2026-11-13. Journal data checked against INSPIRE, arXiv identifiers against arXiv. The comparison with
`notes/warp/phantom_source.md` and `notes/warp/mass_bounds.md`, the search log and what was not read are in
`notes/warp/prior_work.md`. Below, u is the 4-velocity of the matter and n the Eulerian normal.

### Matter sources
- **Santos-Pereira, Abreu & Ribeiro 2020**, arXiv:2008.06560, EPJC 80, 786. Dust T = μ u u with u = n, for a general
  shift β(t, x, y, z). The field equations force μ = 0 and β independent of y and z, which is vacuum with β obeying the
  inviscid Burgers equation (Sec. IV). Does not treat u ≠ n or any other matter.
- **Santos-Pereira, Abreu & Ribeiro 2021a**, arXiv:2101.11467, EPJC 81, 133. Perfect fluid with u = n. The field
  equations give p = 3μ (Eq. 3.4). With p constant, either β is linear in y with μ ≤ 0, or β is complex, or the
  solution is vacuum. A "parametrized perfect fluid" with extra stresses and heat flux admits μ > 0. Does not treat a
  perfect fluid with u ≠ n, or scalar fields.
- **Santos-Pereira, Abreu & Ribeiro 2021b**, arXiv:2108.10960, EPJ Plus 136, 902. Perfect fluid with Λ and u = n
  (Eqs. 3.2–3.9). Every solution is vacuum (Λ = μ = p = 0, Burgers) or has β depending on (t, y) or (t, z) only and
  complex. Same group, read only in abstract and setup: charged dust with u = n (arXiv:2102.05119, GRG 53, 23), the MG16
  summary (arXiv:2111.01298), Burgers reductions (arXiv:2510.11836, EPJ Plus 141, 921), Darmois matching to Minkowski
  (arXiv:2512.12541, EPJC 86, 46) and the thesis that collects them (arXiv:2508.20348).
- **Abellán, Bolívar & Vasilev 2023a**, arXiv:2302.13826, GRG 55, 60. Anisotropic fluid with u = n and principal axes
  x, y, z, general β. The field equations give p_x = 3ρ = ±(3/2)(p_y − p_z) (Eqs. 30–31). Dust, an isotropic fluid, or
  equal tangential pressures force ρ = p_x = 0 and β independent of y or of z, for Cartesian and cylindrical warps.
  Does not treat u ≠ n. Companion: spherical warp metrics with fluids, numerical (arXiv:2305.03736, EPJC 83, 7;
  abstract and setup read).
- **Béatrix-Drouhet 2020**, arXiv:2012.09941 (preprint, "submitted to JBIS"). Fits the Alcubierre G_ab (σ = 8, R = 1)
  with mixtures of two or three anisotropic fluids in boosted and rotated tetrads, by least squares on a grid. Lists
  seven nonzero cylindrical components (tt, XX, tX, ρρ, ρX, ρt, ϕϕ), no ϕ cross terms. Says one fluid (four unknowns,
  seven equations) gives an overdetermined system that is "usually inconsistent", without a proof.
- **Clough, Dietrich & Khan 2024** [CDK24], arXiv:2406.02466, Open J. Astrophys. 7. Initial matter is the whole
  G_ab/8π, with anisotropic stress and momentum, relaxed in time towards a stiff fluid at rest in the Eulerian frame. A
  fluid with a general rest frame is left for later work (Sec. 2.2, Conclusions).
- **Lentz 2021**, arXiv:2006.07125, CQG 38, 075015. Claims a positive-energy soliton sourced by an Einstein–Maxwell
  plasma. Disputed by Santiago, Schuster & Visser 2022 (App. C), by Celmaster & Rubin 2025 (arXiv:2511.18251, preprint:
  the Eulerian energy density is negative in regions) and by Barzegar, Buchert & Vigneron 2026 (Errors 16–17). Abstracts
  and the critiques read, the claim not checked here.
- **Fell & Heisenberg 2021**, arXiv:2104.06488, CQG 38, 155020. Zero-vorticity drive with non-negative Eulerian energy;
  Santiago et al. 2022 show other observers still see the WEC fail. Abstract read.
- **Buchert & Frackowiak 2026**, arXiv:2605.03653, Universe 12, 132. Synge's G-method for Natário-class metrics;
  proposes tilted fluid flows for later work. Abstract read.
- **González-Díaz 2007**, Phys. Lett. B 657, 15, doi:10.1016/j.physletb.2007.09.041. Phantom energy accreted onto a
  two-dimensional superluminal bubble, a test fluid, not a source of the metric. *Abstract only (INSPIRE), no arXiv
  version.*
- Modified sources, abstracts read: DeBenedictis & Ilijić 2018, arXiv:1807.09745, CQG 35, 215001 (Einstein–Cartan
  spin can restore the WEC and NEC); Huey 2024, arXiv:2311.07193, CQG 41, 135007 (thin membranes with perfect-fluid-like
  surface stress, not compact).

### Energy conditions
- **Pfenning & Ford 1997**, arXiv:gr-qc/9702026, CQG 14, 1743. Total Eulerian energy E = −(v²/12)∫ r² f′² dr (Eq. 26),
  with a thin-wall estimate (Eq. 28), and quantum-inequality bounds. No momentum density.
- **Natário 2002**, arXiv:gr-qc/0110086, CQG 19, 1157. Theorem 1.7: a non-flat warp drive (unit lapse, flat slices)
  violates the WEC or the SEC. Introduces the zero-expansion drive.
- **Lobo & Visser 2004a**, arXiv:gr-qc/0406083, CQG 21, 5871, and **2004b**, arXiv:gr-qc/0412065 (Spanish Relativity
  Meeting proceedings). For null k along ± the direction of motion, the exact T_ab k^a k^b is −(v²/8π)|∇⊥f|² ±
  (v/8π)∇⊥²f (2004a Eq. 15), so the NEC fails at arbitrarily small v. The ± average is 4ρ (Eq. 16), which is T_zz = 3ρ
  in the Eulerian frame. Transverse null directions likewise (Eqs. 18–19). In linearized theory T_0i = O(v) and
  T_ij = O(v²) (2004a Sec. III), and boosted observers see ρ < 0 at O(v) (2004b Eq. 8). The O(v) term integrates to
  zero, so they quantify the net violation as O(v²) ≈ −v²R²σ. They state, without computing it, that the ADM mass is
  zero. They do not integrate |J| or split T.
- **Alcubierre & Lobo 2017**, "Warp drive basics", arXiv:2103.05610, Fundam. Theor. Phys. 189, 257–279. Review; repeats
  the NEC computation of Lobo & Visser (Eqs. 12–13).
- **Santiago, Schuster & Visser 2022**, arXiv:2105.03079, PRD 105, 064038. For unit lapse and flat slices:
  ρ = (K² − tr K²)/16π and the Eulerian flux f_i = (1/16π)(∇×∇×v)_i with v = −β (Eqs. 4.7–4.8), which is the J of
  mass_bounds.md item 1. All such drives violate the WEC, SEC and DEC, and the NEC under subsidiary conditions. For
  Alcubierre they cite G_zz = 3G_nn, "established by explicit computation" with no geometric reason given (Sec. 7.1).
  Zero-vorticity drives are Hawking–Ellis Type I. No algebraic type for Alcubierre, no matter model.
- **Helmerich et al. 2024**, arXiv:2404.03095, CQG 41, 095009 (Warp Factory), and **Fuchs et al. 2024**,
  arXiv:2405.02709, CQG 41, 095013. Numerical NEC, WEC, SEC and DEC over sampled observers. Fuchs et al. add a shift to
  a matter shell of positive ADM mass and report no violations; a larger shift adds momentum flux, which the energy
  density must exceed, and the bound is left open (Sec. 5.2). Grep-level reading. Barzegar et al. 2026 (Error 18) say
  the shell metric does not solve the Einstein equations.
- **Rodal 2026**, arXiv:2512.18008, GRG 58, 1. Eigenvalues of G^a_b for one parameter set. Alcubierre is Hawking–Ellis
  Type IV along the whole axis of motion and over the equatorial cut, except an annulus of width ≈ 2Δ about the wall,
  which is Type I (Sec. 4.2). Type IV implies NEC violation. Builds an irrotational drive that is Type I everywhere.
  Does not ask which matter could be the source.
- **Le 2026**, arXiv:2602.18023, CQG 43, 215013. Observer-complete NEC, WEC, SEC and DEC tests as matrix inequalities.
  For flat-slice, unit-lapse drives with β ∝ v, ρ_n and S_ij ∝ v² and |j| ∝ v exactly. When ĵ is a principal axis of
  S_ij, the (n, ĵ) block has complex eigenvalues iff 2|j| > |ρ_n + S_∥|, and the null witness k = n ± ĵ gives
  ρ_n + S_∥ − 2|j|, linear in |j| (Eqs. 18–19). Numerically the Alcubierre wall (f ∈ [0.1, 0.9]) has a Type I
  fraction of 1.3 % at v = 0.5, rising to 24.9 % at v = 2.5, the rest Type IV (Sec. 3, Tables 2–3). The integrated
  negative Eulerian energy scales exactly as v².
- **Barzegar, Buchert & Vigneron 2026**, arXiv:2602.16495 (preprint). Classification and errata of warp-drive papers.
  Corollary IV.37: asymptotically flat flat-slice drives violate the DEC, by the positive energy theorem. See also the
  ADM entry below.
- **Bolívar, Vasilev & Abellán 2026**, arXiv:2609.36211 (preprint). Sharp lower bound on the Lobo–Visser Eulerian
  negative energy of Natário drives, ∝ v²R⁴/Δ³. Eulerian energy only. Abstract read.
- **Bobrick & Martire 2021**, arXiv:2102.06824, CQG 38, 105009. Eulerian energy and momentum integrals; every warp drive
  is a shell of matter; a positive-energy subluminal class. Abstract read, the rest at grep level.

### ADM quantities and the Witten spinor
- **Schuster, Santiago & Visser 2023**, arXiv:2205.15950, GRG 55, 14. States that the Alcubierre and zero-expansion
  drives have identically zero ADM mass, and builds a zero-vorticity drive with ADM mass M (flow ~ √(2M/r)). ADM momentum
  not discussed.
- **Barzegar, Buchert & Vigneron 2026**, arXiv:2602.16495 (preprint). Theorem IV.19: flat-slice ("R-Warp") drives that
  are asymptotically flat have E_ADM = 0, P_ADM ≠ 0 in general, and p_ADM = 0 for Alcubierre. Theorem IV.20: with the
  DEC such a drive is Minkowski. No spinor argument.
- **Frauendiener, Nester & Szabados 2011**, arXiv:1105.5008, CQG 28, 185004. The 2-spinor Witten operator is the
  intrinsic Dirac operator plus ½(tr K)t (Eq. 2.1), so only tr K enters. On maximal, intrinsically flat slices the
  constant spinors solve the Witten equation (conformally constant ones on maximal conformally flat slices). Studies
  boundary conditions for nowhere-vanishing solutions. Nothing on tr K ≠ 0.
- **Frank & Loss 2022**, arXiv:2201.03610 (JEMS 28, 305 (2026) *(unverified)*). Zero modes of γ·(−i∇ − A)ψ = 0 on ℝ^d.
  They recall the standard argument (Sobolev, Hölder, diamagnetic inequality) that ψ ≢ 0 needs ‖A‖²_{L^d} ≥ S_d, and
  prove the sharp ‖A‖²_{L^d} ≥ d S_d/(d − 2), which is 3S₃ in three dimensions (Theorem 1). Remark 4 extends it to
  |γ·∇ψ| ≤ |A||ψ| with ψ ∈ L^{2d/(d−1)}. S₃ = 3(π/2)^{4/3} is the constant of mass_bounds.md item 5. Not about gravity.
  Earlier: Frank & Loss, arXiv:2012.13646 (magnetic-field version, abstract read).

## Cosmology / observation
- Farnes 2018, A&A 620, A92; critique Socas-Navarro 2019, A&A 626, A5 (large-scale runaway).
- Trivedi & Loeb 2026, arXiv:2605.10976 — GW signatures of negative-mass binaries.
- ALPHA-g 2023: antihydrogen falls down, g_H̄ ≈ 0.75 ± 0.29 g; ALPHA-g, AEgIS, GBAR target 1%.

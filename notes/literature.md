# Literature notes

## Foundations
- **Bondi 1957**, Rev. Mod. Phys. 29, 423 — inertial / passive / active mass; runaway pair from Weyl metrics
  (negative-mass source cancels the conical defect). *Read from a PDF copy for the novelty check
  ([novelty/RESULTS.md](novelty/RESULTS.md) §C).* Chapel Hill version:
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
- **Gleiser & Dotti 2006**, gr-qc/0604021, Class. Quantum Grav. 23 (2006) 5063–5078, doi:10.1088/0264-9381/23/15/021
  (arXiv listing and INSPIRE) — kinematic singularity at r = 3μ/(2λ); exact smooth unstable modes. Two arXiv versions,
  read 2026-11-23 (#149). v1 (5 Apr 2006) has no Lemma 1, (nv1), (im) or (cz); its expansion (exp) of the ζ-ODE at r_s
  forces ζ′(r_s) = −2λ²/(3μ(2λ+3)) ζ(r_s). v2 (29 May 2006) prints Lemma 1's (cz) (eq. 16) with +2λ²/(3μ(2λ+3)), the
  opposite of its own (exp) and (RW17). Its (nv1) (eq. 14) has the K that follows from (3)–(4), an H₁ 2/r times the
  derived one, and the sign of H₂'s ζ term reversed. The paper's explicit mode (RW24a) fits the derived map, not the
  printed H₁ and H₂ (#146, #149; phase2/dg_proof/zerilli_gauge_rs.md). Journal text not read (publisher HTTP 403),
  *(unverified)*.
- **Dotti & Gleiser 2008**, arXiv:0809.3615 — intertwiner Ψ̂ = (∂ₓ − ψ₀'/ψ₀)Ψ_z; unique self-adjoint
  extension; instability excitable by generic data. Corollary: linearly unstable. Its eq. (61) prints
  −ℓ(ℓ+1)/(4|M|x) for the x⁻¹ term of V̂; the exact coefficient is 0 (phase2/dg_proof/RESULTS.md, L2/L3; #142).
  Three arXiv versions (21 Sep 2008, 19 Jun 2009, 17 Sep 2009); Class. Quantum Grav. 26 (2009) 215002,
  doi:10.1088/0264-9381/26/21/215002 (listing and INSPIRE; journal text not read, HTTP 403, *(unverified)*). Lemma 1
  (v2, v3; unnumbered (crc) in v1): a metric perturbation is smooth iff Ψ_z is C^∞ off r_s = −3M/λ with a Laurent
  expansion Σ_{j≥−1} c_j(r − r_s)^j, c₀ = λ²c₋₁/(3M(3+2λ)). With ζ = (2λr − 3μ)Ψ_z, μ = −2M, this is
  ζ′(r_s) = −2λ²/(3μ(2λ+3)) ζ(r_s), the sign #146 derives. v1 says (crc) "fixes a sign error in the analogous equation
  (16)" of DG06 (#149).
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
- **Torii & Shinkai 2013** (TS13), "Wormholes in higher dimensional space-time: Exact solutions and their linear stability
  analysis", Phys. Rev. D 88, 064027, arXiv:1309.2058, doi:10.1103/PhysRevD.88.064027. Read 2026-11-26 (#202): arXiv v1
  (9 Sep 2013, 6 pages, the only arXiv version) as page image, text layer and TeX source, which agree on Table I. Reference
  checked against arXiv, INSPIRE (recid 1253353) and Crossref (volume 88, issue 6, article 064027, issued 2013-09-13).
  Journal text not read (publisher HTTP 403): whether its Table I differs from v1 is *(unverified)*.
  - **Contains.** The ghost-scalar (ε = −1, massless) wormhole in n = D spacetime dimensions, static, spherical and
    reflection symmetric: f ≡ 1, R′ = √(1 − (a/R)^{2(n−3)}), Ellis's R = √(r² + a²) at n = 4, throat radius scaled to a = 1.
    Linear stability of the spherical modes e^{iωt}: a gauge-invariant ψ₁ with −ψ₁″ + Vψ₁ = ω²ψ₁ (V singular at the
    throat), regularized with the zero mode ψ̄₁ = 1/(R₀^{(n−4)/2}R₀′) into −Ψ₁″ + W(r)Ψ₁ = ω²Ψ₁ (3.17), W(r) =
    −[3(n−2)²R₀^{−2(n−3)} − (n−4)(n−6)]/(4R₀²) (3.18). Table I lists the one negative ω² found at n = 4 … 11, 20, 50, 100,
    each printed to 15 significant digits, and says the time scale of the instability shortens as n grows.
  - **Does not contain.** How the eigenvalues were computed (no scheme, grid, domain or tolerance; the acknowledgments
    name only the computers). Any uncertainty or statement of precision for Table I. A proof that there is exactly one
    negative eigenvalue (the abstract says "all higher-dimensional wormholes"; Table I has 11 values of n). Values at other n,
    any formula or large-n behavior, tabulated eigenfunctions (Fig. 3 only), nonlinear evolution (deferred to later papers).
  - **Table I against #195's crate** (`notes/wormhole_largeD/ts13_digits.py`, output `ts13_digits_output.txt`). W²a² = −ω²
    (a = 1, e^{Wt} = e^{iωt}). The crate's quad-double values are rounded to 15 digits here.

    | n | ω² as printed | crate W²a² | crate − (−ω²) |
    |---|---|---|---|
    | 4 | −1.39705243371511 | 1.39705243433351 | +6.18e−10 |
    | 5 | −2.98495893027790 | 2.98495893035176 | +7.39e−11 |
    | 6 | −4.68662054299460 | 4.68662054311946 | +1.25e−10 |
    | 7 | −6.46258414126318 | 6.46258414136979 | +1.07e−10 |
    | 8 | −8.28975936306259 | 8.28975936307248 | +9.89e−12 |
    | 9 | −10.1535530451867 | 10.1535530449455 | −2.41e−10 |
    | 10 | −12.0442650147438 | 12.0442650144979 | −2.46e−10 |
    | 11 | −13.9552091676647 | 13.9552091674964 | −1.68e−10 |
    | 20 | −31.5751101285105 | 31.5751101283585 | −1.52e−10 |
    | 50 | −91.3457759137153 | 91.3457759136296 | −8.57e−11 |
    | 100 | −191.283017729717 | 191.283017730359 | +6.42e−10 |

    The printed values agree with the crate to 9.35 … 12.03 significant digits, at no row to the last printed digit.
    `one_mode.log` (N4: shooting for TS13's own (3.18) potential, which `ts13_digits.py` checks equals its U exactly; no code
    shared with the crate) equals the crate at all 12 digits it prints at n = 4, 10, 20, 50. Table I minus that value is
    −6.2e−10, +2.4e−10, +1.1e−10, +1.2e−10 there. That is **numerical evidence** that Table I's digits beyond about the 9th
    to 12th significant are not accurate at those four n. At the other seven rows the evidence is the crate alone (its qd and
    od lanes agree to 60 digits but share code). growth.csv's D = 4 … 11 rows are Table I rounded to 12 decimals, so they
    carry that error; no csv digit there goes beyond the printed ones. README.md's quoted [TS13] values are roundings of
    Table I.

## Warp drives: matter sources, energy conditions, ADM quantities
Read for #60 on 2026-11-13. Journal data checked against INSPIRE, arXiv identifiers against arXiv. The comparison with
`notes/warp/phantom_source.md` and `notes/warp/mass_bounds.md`, the search log and what was not read are in
`notes/warp/prior_work.md`. Below, u is the 4-velocity of the matter and n the Eulerian normal.

### The metric
- **Alcubierre 1994**, arXiv:gr-qc/0009013, Class. Quantum Grav. 11, L73–L77, doi:10.1088/0264-9381/11/5/001 (checked on
  arXiv). Defines ds² = −dt² + (dx − v_s f(r_s) dt)² + dy² + dz² with the tanh profile, and gives the Eulerian energy
  density −(v_s²/32π)((y² + z²)/r_s²) f′² < 0 (Eq. 19; the arXiv text prints r_c where r_s is meant, as #57 found).
  Concludes that exotic matter is needed. Names no matter model. Reproduced in `notes/warp/alcubierre_bounds.py`.

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
  violates the WEC or the SEC. Introduces the zero-expansion drive. Secs. 1 and 2 read in the arXiv v3 text for #66,
  whose closed-form μ equals Sec. 2's ρ ([warp/natario.md](warp/natario.md)).
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

## Universe from collapse: ghost wormholes, phantom black universes, black bounces, black-hole interiors
Read 2026-11-27 for #211, from local copies of the arXiv texts. What was read of each paper, its version, and the
verification of every reference against INSPIRE and arXiv are in negative-mass-raw `I211-rpi5/README.md`. Per entry: the
claim in the paper's words; the mechanism, with the energy condition it violates and the matter; what the paper contains
(a solution, its formation, its stability); what it does not contain. "Abstract only" means the full text was not read.
The ranked list of claims that a computation here could test is [baby_universes.md](baby_universes.md). The ghost
wormhole's thresholds and exponents are in `critical_ghost.md` §1.2–1.3 and are not repeated; [GGS09a] and [TS13] are in
its reference list.

### The ghost wormhole's expanding branch
- **[SH02] Shinkai & Hayward 2002**, "Fate of the first traversible wormhole: black-hole collapse or inflationary
  expansion", PRD 66, 044005, arXiv:gr-qc/0205041 v2, doi:10.1103/PhysRevD.66.044005 (INSPIRE and arXiv checked). Read
  in full.
  - Claim: "The wormhole throat suffers a bifurcation of horizons and either explodes to form an inflationary universe
    or collapses to a black hole, if the total input energy is respectively negative or positive." "The inflationary
    case provides a mechanism for inflating, to macroscopic size, a Planck-sized wormhole formed in space-time foam." In
    Sec. III B: "This exponential expansion, combined with the horizon structure, indicates that the wormhole has
    exploded to an inflationary universe. The two trapping horizons become the cosmological horizons." In Sec. IV: "for
    one of the simplest classical exotic matter fields, wormholes can naturally inflate. The problem is only to stop the
    inflation proceeding indefinitely."
  - Mechanism: Einstein gravity with a massless ghost Klein–Gordon field, R_ab = 2∇ψ∇ψ − 2∇φ∇φ (Eq. 2.1), spherical
    symmetry, dual-null evolution that covers both universes. For null k, R_kk = −2(k·∇φ)² ≤ 0 from (2.1) (by hand): the
    NEC fails wherever the ghost varies along k. The paper's words: the Bondi energy-loss property that holds for NEC
    matter "is reversed for the ghost field".
  - Contains: evolution to x⁺ = 20 of the Ellis wormhole with Gaussian pulses in the ghost field or in a normal field.
    The sign of the initial Bondi energy E₀ decides black hole (E₀ > 0) or expansion (E₀ < 0), with the final structure
    "judged at x⁻ = 10 by the horizon structure". For the expanding case a fit r/a = 1 + b₄ exp(H(τ − b₅)) of the areal
    radius of the diagonal x⁺ = x⁻ against its proper time, with Ha = 1.56, 1.19, 1.13, 1.12 for pulse amplitudes
    10⁻¹…10⁻⁴ (Table I) and H ∼ 1.1/a (3.7). The horizon bifurcation (Figs 4–5). The paper assumes the static wormhole
    is linearly stable ("Accepting that the static wormhole is stable to linear perturbations [11, 12], we seem to have
    discovered a non-linear instability"); [GGS09a] found the linear mode.
  - Does not contain: the matter content of the "universe" (no ρ, p or w of the expanding region, no comparison with a
    de Sitter or Kantowski–Sachs solution); a run past x⁻ = 10; a mechanism or condition that ends the expansion; the
    energy bookkeeping of the expanding region. The data are asymptotically flat at both ends and the static massless
    wormhole has zero Bondi energy; that what "inflates" is the region round the throat is an inference from the setup,
    not a statement of the paper. The "throat" is the diagonal, which the paper says "loses its specific geometrical
    meaning" once perturbed (note 17). How a Planck-sized wormhole forms in foam is not treated.
  - Reproduced here (#39, `wormhole_evolution/README.md`): both fates and the sign rule of E₀, and the linear rate to
    3e-5. Its expansion probe saw no coordinate limit to t = 42 and H on the diagonal rising 0.75 → 2.9, not interpreted
    (#213).
- **[GGS09b] González, Guzmán & Sarbach 2009**, "Instability of wormholes supported by a ghost scalar field. II.
  Nonlinear evolution", CQG 26, 015011, arXiv:0806.1370 v2, doi:10.1088/0264-9381/26/1/015011 (INSPIRE and arXiv
  checked). Read: abstract, §I, §V–VII; §II–IV skimmed.
  - Claim: "In the expanding case, the wormholes seem to undergo an exponential expansion, at least during the run time
    of our simulations." "In the expanding case, we do not observe any apparent horizons and our simulations suggest
    that the wormhole expands forever."
  - Mechanism: the same ghost field in Cauchy form, throat radius b, the γ₁ family (opposite-sign ADM masses for γ₁ ≠
    0), two gauges (λ = 0, 1). The NEC fails as in [SH02].
  - Contains: for reflection-symmetric data the areal radius at x = 0 fits r(τ) = 1 + e^{(τ−p₁)/p₂}, p₂ → 0.85 as ε_c →
    0 "in good agreement with the prediction from perturbation theory" (Fig. 10, proper-time axis to 8). No apparent
    horizon in any expanding run. The scalar "does not seem to dissipate away but actually gains in strength near the
    point of reflection symmetry". The throat is slicing-dependent: in the λ = 1 gauge r shows two minima, "a pure gauge
    effect". The amplitude scan with the second collapse window is in `critical_ghost.md` §1.2.
  - Does not contain: the word "inflation" in its text (it appears only in the title of [SH02] in the reference list); a
    Hubble rate, an equation of state or any cosmological reading of the expansion; the energy of the expanding state; a
    run beyond the proper time of Fig. 10.

### Phantom black universes
- **[BF06] Bronnikov & Fabris 2006**, "Regular phantom black holes", PRL 96, 251101, arXiv:gr-qc/0511109 v1,
  doi:10.1103/PhysRevLett.96.251101 (INSPIRE and arXiv checked). Read in full, PDF and TeX.
  - Claim: "The Penrose diagram of a regular BH is Schwarzschild-like, but the singularity at r = 0 is replaced by a de
    Sitter infinity, which gives a hypothetic BH explorer a chance to survive. Such solutions also lead to the idea that
    our Universe could be created from a phantom-dominated collapse in another universe, with KS expansion and
    isotropization after crossing the horizon." In the body: "One may also speculate that our Universe could appear from
    collapse to a phantom BH in another, 'mother' universe and undergo isotropization (e.g., due to particle creation)
    soon after crossing the horizon."
  - Mechanism: a minimally coupled scalar with kinetic sign ε = −1 and a potential V(φ), static spherical metric ds² = A
    dt² − dρ²/A − r²dΩ². "Since it does not respect the usual energy conditions": Eq. (5), 2r″/r = −εφ′², gives r″ ≥ 0,
    so a minimal sphere is allowed, which ε = +1 forbids (no-go A). Eq. (9)–(14) is the inverse-problem example r = √(ρ²
    + b²), φ = √2 arctan(ρ/b), two constants c and ρ₀. It is asymptotically flat at ρ → +∞ iff 2bc = −πρ₀, m = ρ₀/3. For
    m > 0 it is a regular black hole with a de Sitter far end beyond the horizon. The paper says every model of its
    classes needs "a potential with at least two zero-slope points".
  - Contains: a classification of 16 classes of regular static solutions of a phantom scalar; the explicit family; the
    claim that de Sitter is a global attractor of isotropic phantom cosmologies with V bounded above (Faraoni), with
    "Very probably this conclusion extends to KS cosmologies after isotropization". The universe-from-collapse idea is a
    stated speculation.
  - Does not contain: a collapse or any dynamics, so no formation; a perturbation or stability analysis (its
    introduction says that "in a classical setting, a massless phantom field even shows a more stable behavior than its
    usual counterpart [12, 13]", and [12, 13] are the works [SH02] cites for linear stability of the massless wormhole,
    where [GGS09a] later found a growing mode); isotropization, which is named as possible "due to particle creation"
    and not computed.
  - Printed equations, from #212 (`black_universe/README.md` §1): (8) prints A′(h) = (ρ₀ − h)/r²(h), while (7), B′ ≡
    (A/r²)′ = 2(ρ₀ − h)/r⁴ with A(h) = 0, gives A′(h) = 2(ρ₀ − h)/r²(h). Checked here by hand from the printed (7); #212
    S3 confirms it to 40 digits. The printed (8) halves the surface gravity κ = |A′(h)|/2. #212 also finds the closed
    forms (9)–(14) satisfy (3)–(6) with residual 0 (S1, S2, *checked symbolically*).
- **[BMD07] Bronnikov, Dehnen & Melnikov 2007** (author order corrected 2026-11-27 against arXiv gr-qc/0611022 and Crossref), "Regular black holes and black universes", GRG 39, 973,
  arXiv:gr-qc/0611022 v2, doi:10.1007/s10714-007-0430-6 (INSPIRE and arXiv checked). Read: abstract, §1, concluding
  paragraphs.
  - Claim: "It also looks possible that our Universe has originated from phantom-dominated collapse in another universe
    and underwent isotropization after crossing the horizon." "We are thus facing one more mechanism of universes
    multiplication, in addition to the well-known mechanism existing in the chaotic inflation scenario."
  - Mechanism: static distributions of phantom matter in scalar-tensor gravity and k-essence, with "necessary conditions
    for the existence of black-universe solutions". The NEC fails as in [BF06].
  - Contains: a comparison of types of regular static black holes; the name "black universes" for those with an
    expanding Kantowski–Sachs region beyond the horizon. Does not contain: formation, stability or the isotropization
    computed.
- **[BKZ12] Bronnikov, Konoplya & Zhidenko 2012**, "Instabilities of wormholes and regular black holes supported by a
  phantom scalar field", PRD 86, 024028, arXiv:1205.2224 v3, doi:10.1103/PhysRevD.86.024028 (INSPIRE and arXiv checked).
  Read: abstract, §I–II, §V C, §VI; printed (36)–(41) from the TeX; §III (axial formalism) and the quasinormal-mode
  computation not read.
  - Claim: "we have shown that all configurations under study are unstable under spherically symmetric perturbations,
    except for a special class of black universes where the event horizon coincides with the minimum of the area
    function." Conclusion 2: black universes "are shown to be stable only in the special case where the horizon
    coincides with the minimum of the area function R(r) (the parameters: c = −1, b = 3πm/2) and unstable for c ≠ −1."
    On approaching c = −1: "we do not observe a decrease in the growth rate as expected ... This can be an indication
    that the only parameter for which a black universe can be stable is c = −1."
  - Mechanism: the [BF06] family with a phantom scalar, linearized about the static solution: general axial
    perturbations and the monopole polar mode, a Schrödinger problem whose potential has an infinite wall at a throat,
    regularized numerically, plus time-domain integration.
  - Contains: M-AdS wormholes unstable (regularized potentials, Figs 1–2); black universes unstable except c = −1, by
    potentials and time-domain runs at several c; axial modes stable in the static region; axial quasinormal frequencies
    at c = −1 (Table I).
  - Does not contain: nonlinear evolution (listed as open: "a full nonlinear analysis of perturbations in all relevant
    cases"); formation; polar ℓ ≥ 1; a count of growing modes or any growth rate (no rates are given); the
    Kantowski–Sachs region.
  - #212 result (`black_universe/README.md`, *numerical evidence*): the verdict agrees at all 20 members tested.
    Nineteen members (c = −10…−0.01) have exactly one growing radial mode, with rates 0.09…1.17 per b; c = −1 has none.
    The "class" is one member of a one-parameter family, and the rate tends to κ = 2/π per b as c → −1, not to zero.
  - Printed equations, from #212 (§2): the derived wave equation differs from the printed (36)–(38) and agrees with
    (40), (41). Printed (36): time-derivative coefficient e^{−2μ₂−2ν} (= 1 for e^{2ν} = A, e^{2μ₂} = A⁻¹; the paper's
    definitions give 1/A², by hand here) and δφ′ coefficient ν′ − 2μ₃′ + μ₂′ (derived: ν′ − μ₂′ + 2μ₃′ = A′/A + 2R′/R).
    Printed (37): first term ε(e^{2μ₃} − V)φ′²/μ₃′² (derived: ε(V − e^{−2μ₃})φ′²/μ₃′²) and V_φ term −2φ′V_φ/μ₃′
    (derived: +2φ′V_φ/μ₃′). Printed (38): δφ = Ψe^{μ₃} (derived: Ψ/R). Only the first row was rechecked here; the others
    are #212's sympy results.
- **[BBS12] Bolokhov, Bronnikov & Skvortsova 2012**, "Magnetic black universes and wormholes with a phantom scalar", CQG
  29, 245006, arXiv:1208.4619 v1, doi:10.1088/0264-9381/29/24/245006 (INSPIRE and arXiv checked). Read: abstract,
  opening of §1, the stability paragraph of the conclusions; the rest by search.
  - Claim: "The black universe models make possible a cosmological scenario where a phantom-dominated gravitational
    collapse in some 'mother' universe creates our Universe whose expansion begins from a horizon, and the next stages
    are isotropization and de Sitter inflationary expansion."
  - Mechanism: a phantom scalar with a potential plus a magnetic field. Contains: explicit globally regular solutions
    with 0–4 horizons. Does not contain: formation or a stability computation. It says that "Most of the known
    scalar-vacuum wormhole and black-universe solutions proved to be unstable under radial perturbations ... and it is
    of interest to find out whether or not they can be stabilized by electric or magnetic fields", left for later.
- **[BCK18] Blázquez-Salcedo, Chew & Kunz 2018**, "Scalar and axial quasinormal modes of massive static phantom wormholes",
  PRD 98, 044035, arXiv:1806.03282, doi:10.1103/PhysRevD.98.044035 (INSPIRE and arXiv checked). Read: the TeX source, §V
  (radial, scalar and axial results) and its figures' data (the EPS files).
  - Claim: "scalar, axial and radial perturbations of massive static spherically symmetric wormholes supported by a phantom
    scalar field", quasinormal modes by direct integration and WKB. The unstable radial mode at M = 0 has ω_I = −1.182 in
    its convention, i.e. W = 1.182 per r0 (#212: 1.18197); the scalar and axial modes are damped.
  - Contains: figures only, no tables. From the EPS data at M = 0 (r0 = 1; read from the plotted points, three digits): the
    scalar branches ω_R = 0.684 (ℓ = 0), 1.573 (ℓ = 1), 2.547 (ℓ = 2), and ω_I(ℓ = 1) = 0.530 (their sign convention, damped).
    #218 reproduces the ℓ = 1 value as the quasinormal mode of the *test*-scalar potential 2/R² + 1/R⁴: 1.572708 − 0.529701i.
    Does not contain: polar gravitational perturbations.
- **[ABCKY23] Azad, Blázquez-Salcedo, Chew, Kunz & Yeom 2023**, "Polar modes and isospectrality of Ellis–Bronnikov wormholes",
  PRD 107, 084024, arXiv:2212.12601, doi:10.1103/PhysRevD.107.084024 (INSPIRE and arXiv checked). Read: the TeX source and the
  PDF text of the arXiv version of 27 Dec 2022 (equation numbers below are that version's), all sections and Table 1.
  - Claim: "For massless wormholes the perturbations decouple, and we obtain two identical master equations for the scalar
    and gravitational modes, which moreover agree with the master equation for the axial modes. Consequently there is
    isospectrality with threefold degenerate modes." For C = 0, ℓ ≥ 2: Z″ + (ω² − V)Z = 0, V = ℓ(ℓ+1)/(r² + r0²) − 3r0²/(r² + r0²)²
    ((52), (53)). Table 1, row M/r0 = 0: ω r0 = 1.74 − 0.305i, 2.95 − 0.410i, 4.08 − 0.449i for ℓ = 2, 3, 4.
  - Contains: the polar system in Regge–Wheeler gauge for C ≠ 0 and C = 0, ℓ ≥ 2 numerically; ℓ = 1 and ℓ = 0 by a scalar
    equation, (69) and (72), with no numbers. Does not contain: a count of growing polar modes, members with a potential
    (black universes, M-AdS), or any value for ℓ = 1.
  - #218 (`black_universe/polar.md`, *checked symbolically* and *numerical evidence*): the derived ℓ ≥ 2 system reproduces (44)–(46) and
    the master equation (52), (53), and Table 1's M = 0 row (1.737692 − 0.305138i, 2.952415 − 0.409984i, 4.076253 − 0.449105i).
    Three disagreements with the printed text, reported and not corrected (polar.md Sec. 5): (47) has the sign of its
    (ℓ+2)(ℓ−1)K term opposite to the derivation; (69) has the sign of its 2/R² term opposite to the test scalar's, and
    neither is the derived ℓ = 1 equation (the gauge-invariant scalar obeys V = 1/R⁴, the Darboux partner of V at L = 2);
    the statement that ℓ = 1 "corresponds to the scalar branch" of [BCK18] holds, in their wording, only for the backreacted
    potential 1/R⁴, not for the test field's 2/R² + 1/R⁴.

### Black bounces
- **[SV19] Simpson & Visser 2019**, "Black-bounce to traversable wormhole", JCAP 02 (2019) 042, arXiv:1812.07114 v3,
  doi:10.1088/1475-7516/2019/02/042 (INSPIRE and arXiv checked). Read: abstract, §1, §2 text on the bounce, §5 opening
  and closing, §8 opening, §9.
  - Claim: the metric ds² = −(1 − 2m/√(r² + a²))dt² + dr²/(1 − 2m/√(r² + a²)) + (r² + a²)dΩ² "neatly interpolates
    between the standard Schwarzschild black hole and the Morris–Thorne traversable wormhole; at intermediate stages
    passing through a black-bounce (into a future incarnation of the universe), an extremal null-bounce ... and a
    traversable wormhole." In §9 the singularity "is converted into a spacelike hypersurface representing a 'bounce'
    into a future incarnation of the universe".
  - Mechanism: a metric ansatz, with no matter model. From the Einstein equations (§5) the stress-energy "clearly
    violates all the classical energy conditions": the NEC fails everywhere except on horizons. The region beyond the
    bounce is "a separate copy of our own universe", not a de Sitter region.
  - Contains: Carter–Penrose diagrams, curvature, the stress-energy and energy conditions, surface gravity, ISCO and
    photon sphere, Regge–Wheeler potentials for spin 0 and 1 test fields on the fixed background. Does not contain: a
    matter source, formation by collapse, stability of the geometry under coupled perturbations.
- **[LRSSV21] Lobo, Rodrigues, Silva, Simpson & Visser 2021**, "Novel black-bounce spacetimes: wormholes, regularity,
  energy conditions, and causal structure", PRD 103, 084052, arXiv:2009.12057 v3, doi:10.1103/PhysRevD.103.084052
  (INSPIRE and arXiv checked). Abstract only. Claim: "specific regular black holes where the 'area radius' always
  remains non-zero, thereby leading to a 'throat' that is either timelike ..., spacelike (corresponding to a 'bounce'
  into a future universe), or null". It analyzes regularity, energy conditions and causal structure of metrics from a
  mass function. The abstract names no matter source, formation by collapse or stability analysis; the text was not read
  beyond a search for stability, which finds only a citation of [BKZ12].
- **[BW22] Bronnikov & Walia 2022**, "Field sources for Simpson–Visser spacetimes", PRD 105, 044039, arXiv:2112.13198
  v1, doi:10.1103/PhysRevD.105.044039 (INSPIRE checked; the arXiv page prints "Phys. Rev. D 15", an apparent typo).
  Read: abstract, §I, the NEC statement, concluding remarks.
  - Claim: "any static, spherically symmetric SV metric can be obtained as an exact solution to the Einstein field
    equations sourced by a combination of a minimally coupled phantom scalar field with a nonzero potential V(φ) and a
    magnetic field in the framework of nonlinear electrodynamics". "SV-like geometries cannot be sourced by a scalar
    field or NED taken separately."
  - Mechanism and condition: "The NEC is necessarily violated in a neighborhood of any wormhole throat or any black
    bounce"; with r = √(x² + a²), r″/r = a²/r⁴ > 0, so the NEC fails in the whole space-time except on horizons. The
    phantom scalar carries the violation.
  - Contains: V(φ) and L(F) for the Schwarzschild-like and Reissner–Nordström-like cases; Carter–Penrose diagrams. Says
    of stability: "Another issue of importance is the stability of SV space-times, which problem cannot be formulated
    until we know their precise dynamic context", and the same geometry may have different stability with different
    sources. Does not contain: stability, formation, perturbations.
- **[PPS25] Pal, Pal & Sarkar 2025**, "Dynamical interiors of black-bounce spacetimes", JCAP 01 (2025) 069,
  arXiv:2410.17150 v1, doi:10.1088/1475-7516/2025/01/069 (INSPIRE checked; the arXiv page has no journal reference).
  Read: abstract, §I, conclusions.
  - Claim: a regular FLRW-like interior "can not be matched smoothly with the exterior black-bounce spacetime through a
    timelike hypersurface, as there always exists a thin shell of non-zero energy-momentum tensor at the matching
    hypersurface". For an initially collapsing interior the numerical solution "crosses the regular transition surface
    at l = 0 ... and goes to a separate copy of our universe".
  - Contains: Israel–Darmois matching, a surface energy density and pressure, an evolution equation for the interior
    scale factor, and the formation of an event horizon inside the cloud. Does not contain: a matter model. "In our
    analysis we have not used any detail of the dynamics of the gravitational theory under consideration", and the
    condition on the surface stress is "in a sense, arbitrary", so the bounce behavior may not hold for other
    conditions.
- **[CK21] Chakrabarti & Kar 2021**, "A wormhole geometry from gravitational collapse", PRD 104, 024071,
  arXiv:2106.14761 v1, doi:10.1103/PhysRevD.104.024071 (INSPIRE checked; the arXiv page gives the DOI only). Read:
  abstract, conclusion, passages on the energy conditions.
  - Claim: "gravitational collapse of a NEC (Null Energy Condition) violating spherically symmetric fluid distribution
    can avoid the formation of a zero proper volume singularity and eventually lead to a Lorentzian wormhole geometry".
    Conclusion: "a time-dependent wormhole geometry can arise in gravitational collapse, provided the collapsing fluid
    violates the relevant energy conditions."
  - Mechanism: a time-dependent metric conformal to the static Simpson–Visser form, with an imperfect fluid read off
    from the Einstein tensor ("So far we have not used any theory of gravity as such"). Does not contain: a matter model
    or field equations that produce the fluid, stability, or an evolution from regular data.

### Prior work on the horizon sector of the black universes (#260)
Read 2026-11-27 for #260, from local copies of the arXiv texts, against four statements of this project: S1 (the horizon
coefficient c₂ of the c = −1 member is a free number of regular Cauchy data, #248), S2 (so the member has a growing
direction of rate κ, against the "stable" verdict), S3 (a non-extremal Killing horizon at an extremum of the areal
radius makes c₂ gauge invariant, #261) and S4 (the closed form `Y = R e^{κx}` of the growing mode, #262). Verdicts and
the comparison of boundary conditions are in `black_universe/horizon_literature.md`; the queries (with hit counts), the
screening of 424 full texts, the depth of each reading and the verification of each reference are in negative-mass-raw
`I260-mac`. Per entry: the text that bears on S1–S4, with its section or equation; what the paper contains; what it does
not. Entries above this subsection ([BF06], [BKZ12], [BBS12], [SV19], [LRSSV21], [BW22]) are not changed; additions to
them are listed first. The papers' full texts are not committed (third-party content); `I260-mac/papers/downloads.tsv`
has each url and sha256.

**Additions to existing entries** (re-read in full for [BF06], [BKZ12]; passages for the others).
- [BF06], read on every page. After Eq. (14): "the throat ρ = 0 is located in the R region if c > −1 , i.e., if 3πm < 2b
  , at the horizon if 3πm = 2b and in the T region beyond it if 3πm > 2b". For ρ₀ = 0: "at c = −1 they merge into a
  double horizon". Item (C), p. 2: the no-hair theorem holds for ε = −1, V ≤ 0 "and r′ > 0 at and outside the event
  horizon. It will be seen below that there exist BH solutions where the condition r′ > 0 is violated." Eq. (8) prints
  `A′(h) = (ρ₀ − h)/r²(h)`; with Eq. (7), `B′ = 2(ρ₀ − ρ)/r⁴`, A(h) = 0 gives twice that (as #212). No erratum or
  related record on INSPIRE (queries `q16`, `q17`); no later paper of 424 repeats Eq. (8).
- [BKZ12], read on every page (reference list excluded). Gauge, §IV, Eq. (35): "Then we can use the gauge freedom and
  fix" `δψ = δμ₃ ≡ 0 ⟹ δR(r, t) = 0`. Horizon condition, §V.A: "The quasinormal boundary conditions for black hole
  perturbations imply pure incoming waves at the horizon and pure outgoing waves at spatial infinity." (Eq. 43, Ψ ∝
  e^{±iωr*}.) §V.C: "Black universes with c ≤ −1, or equivalently 3πm ≥ 2b, have no throat in the static region", so no
  regularization at c = −1; "As we have shown, black universes at c = −1 (b = 3πm/2) are stable against polar monopole
  perturbations"; "From Fig. 5 one can see that for c = −1 the effective potential has a negative gap, however,
  time-domain integration proves that in this case the black universe is stable." The method is cited, not described
  (§V.B, for the wormholes: "we have used the time-domain integration method proposed by Gundlach et al. [35]"); the
  initial data of the c = −1 run are not given. Fig. 5 caption: "The potentials vanish at the horizon." §V.A: "For
  nonzero scalar potentials V (φ), analytical expressions for static perturbations are unknown." Conclusion 2 as quoted
  above. Not in the paper: a growth rate, a count of modes, the region beyond the horizon, the nonlinear problem.
- [LRSSV21] §IV.A: "If a = 2m, the (maximally extended) spacetime has only extremal horizons, so we have a one-way
  wormhole geometry with an extremal null throat". Its introduction restates [BKZ12]'s verdict as motivation ("except
  for a special class of black universes where the event horizon coincides with the minimum of the area function"); no
  stability computation.
- [BW22]: for 2M = a "the minimum of r coincides with a single horizon" (the "one-way wormhole"); no stability analysis
  (as in the entry above).

**The black universes: stability and its gauge**
- **[BFZ11] Bronnikov, Fabris & Zhidenko 2011**, "On the stability of scalar-vacuum space-times", EPJC 71, 1791,
  arXiv:1109.6576 v2, doi:10.1140/epjc/s10052-011-1791-2 (INSPIRE and arXiv checked). Read: abstract, §1, §2.4–2.6,
  conclusions; the rest by patterns.
  - Claim, §2.5: "The gauge δβ = 0 (the same as δr = 0) is suitable for describing the perturbations at any points
    except those where r′ = 0: these are throats and other critical points of r(u)." "putting δr = 0, we forbid
    perturbations of the throat radius, while there is no physical reason for that". §2.6: "with nonzero δr on the
    throat. So it is necessary to pass on to another gauge"; Eq. (50): δr is finite in the gauge δφ = 0 if φ′ ≠ 0.
  - Contains: the S-deformation regularization of a throat's potential and the instability of the V ≡ 0 solutions. Does
    not contain: a horizon at a critical point of r, free Cauchy data, or the black universes' c = −1 member. Bears on
    S1 (nonzero δr at a critical point of r is admissible) as a statement for throats.
- **[Bro18] Bronnikov 2018**, "Scalar fields as sources for wormholes and regular black holes", Particles 1, 56,
  arXiv:1802.00098 v2, doi:10.3390/particles1010005 (INSPIRE and arXiv checked). Read: §2.2, §4.1–4.3, Fig. 5; the rest
  by patterns.
  - §4.1, after Eq. (70): "The gauge δβ = 0 is convenient for calculations but looks doubtful when applied to
    configurations with throats." "we suppose that the throat radius is not subject to perturbations, whereas its
    changes are in general admissible". Eqs. (64)–(68): δφ = ψe^{−η}, η′ = β′ + h′/2h, U with (V − e^{−2β}) and
    +2φ′V_φ/β′ (the signs of #212's derivation, not of [BKZ12]'s printed (36)–(38)).
  - §4.3: "the perturbation δφ should remain finite in the whole space and also vanish at infinities or horizons". Fig.
    5 caption: "Only one of the models under study is stable (the time-domain plot is asymptotically constant), it is
    the model where the throat coincides with the event horizon." §2.2: "if a minimum of r coincides with a horizon,
    then the same reasoning shows that NEC violation is necessary on either side in its neighborhood".
  - Contains: a review of the above; no new computation of the c = −1 member. Bears on S2 (restates the verdict and the
    horizon requirement on δφ).
- **[BBS20] Bronnikov, Bolokhov & Skvortsova 2020**, "Hybrid metric-Palatini gravity: black holes, wormholes,
  singularities and instabilities", Grav. Cosmol. 26, 212, arXiv:2006.00559 v1, doi:10.1134/S0202289320030044 (INSPIRE
  and arXiv checked). Read: abstract, §5.2–5.3, §6.
  - §5.3, item 2(i): "among the black-universe solutions there is a stable subset, in which the horizon coincides with
    the sphere of minimum radius", with [BKZ12] as the reference. §5.2: at a horizon "we have no reason to require there
    anything more than finiteness of δφ". §6: "a conformally mapped black universe with a horizon at the minimum
    radius".
  - Contains: stability results of conformally mapped solutions by mapping [BKZ12]'s boundary conditions. Does not
    contain: a new computation at c = −1. Bears on S2 (restatement) and on the list of solutions with the minimum of r
    on a horizon.
- **[Bro19] Bronnikov 2019**, "Spherically symmetric black holes and wormholes in hybrid metric-Palatini gravity", Grav.
  Cosmol. 25, 331, arXiv:1908.02012 v2, doi:10.1134/S0202289319040030 (INSPIRE and arXiv checked). Read: abstract, §3.5.
  §3.5: "if we adopt for the horizons the same boundary condition δψ = 0 as is used at infinity (which looks
  reasonable)". Bears on the horizon condition of the school.
- **[BBS21] Bronnikov, Bolokhov & Skvortsova 2021**, "Spherically symmetric space-times in generalized hybrid
  metric-Palatini gravity", Grav. Cosmol. 27, 358, arXiv:2109.01794 v1, doi:10.1134/S0202289321040046 (INSPIRE and arXiv
  checked). Read: abstract, Figs. 4 and 6.
  - Fig. 4 caption (Example 2): "At m≈ 0.212, the horizon occurs at x = 0 and coincides with the minimum of r(x)."; Fig.
    6 (Example 3): m ≈ 0.141. Einstein-frame solutions with a phantom scalar and a potential; whether the horizon is
    simple is not stated in the captions. No stability analysis there.
- **[BD11] Bronnikov & Donskoy 2011**, "Black universes with trapped ghosts", Grav. Cosmol. 17, 176, arXiv:1110.6030 v1,
  doi:10.1134/S0202289311020083 (INSPIRE and arXiv checked). Read: §2 passage. "In the intermediate case where a minimum
  of r(u) coincides with a horizon (A = 0), the condition r′′ > 0 holds in its vicinity". Solutions only; no
  perturbations.
- Other follow-ups from the same group, read by abstract and patterns (G): [BBS12] 1208.4619 (no c = −1 stability
  statement beyond the entry above), 0910.4930 (Bronnikov & Donskoy, brane-world black universes; no perturbations),
  1503.02956 (Bronnikov & Korolyov; the gauge remark "The gauge δβ = 0 is technically the simplest one, but causes
  certain problems when applied to wormholes and other configurations with throats"), 1705.05906 (Bronnikov & Korolyov,
  long throats).
- **Restatements of [BKZ12]'s verdict, no new computation** (each read at the passage): [ES13] Eiroa & Sendra 2013,
  "Regular phantom black hole gravitational lensing", PRD 88, 103007, arXiv:1308.5959 v2, doi:10.1103/PhysRevD.88.103007
  (§6: "when the horizon coincides with the throat, since it was shown in Ref. [21] that the solution is stable");
  [BFR20] Bronnikov, Fabris & Rodrigues 2020, "Instability of some k-essence spacetimes", IJMPD 29, 2050016,
  arXiv:1908.09126 v1, doi:10.1142/S0218271820500169; [KMP23] Karakasis, Mavromatos & Papantonopoulos 2023, "Regular
  compact objects with scalar hair", PRD 108, 024001, arXiv:2305.00058 v3, doi:10.1103/PhysRevD.108.024001; [Sim23]
  Simpson 2023, "Excising Curvature Singularities from General Relativity", arXiv:2304.07383 v1 (INSPIRE lists
  doi:10.26686/wgtn.22338298, no journal); [FDV23] Fabris, Daouda & Velten 2023, "Spherically Symmetric Configurations
  in Unimodular Gravity", Universe 9, 515, arXiv:2310.19037 v1, doi:10.3390/universe9120515 (§VI: "the solutions are
  unstable, except for the regular solution in the very special case the minimum of the areal function coincides with
  the horizon"); [GOPV25] González, Olivares, Papantonopoulos & Vásquez 2026, "Observational tests of regular black
  holes with scalar hair and their stability", PRD 113, 124076, arXiv:2510.18647 v2, doi:10.1103/dd4m-zs39; [BGPV26]
  Bécar, González, Papantonopoulos & Vásquez 2026, "Anomalous decay rate and greybody factors for regular black holes
  with scalar hair", PRD 113, 104020, arXiv:2602.16972 v2, doi:10.1103/gksh-s6pn; [GOPV26] González, Olivares,
  Papantonopoulos & Vásquez 2026, "Timelike geodesics of regular black holes with scalar hair", EPJ Plus 141, 1227,
  arXiv:2604.20661 v2, doi:10.1140/epjp/s13360-026-08428-4 (all INSPIRE and arXiv checked). [BGPV26] and [GOPV26] use
  test fields and geodesics on the family, with the horizon at the throat at the critical scalar charge; [BGPV26]: "This
  occurs precisely when the event horizon coincides with the minimum of the areal radius". With [Bro18] and [BBS20], 11
  papers. None examines `Y ≠ 0` on the horizon.

**Simpson–Visser black bounces and other horizons on a throat**
- **[BIK26] Bronnikov, Ibadov & Kunz 2026**, "Instabilities of Simpson–Visser black bounces sourced by a phantom scalar
  and nonlinear electrodynamics", arXiv:2611.08322 v1 (5 Nov 2026; INSPIRE record 3225530 has no journal reference)
  (INSPIRE and arXiv checked). Read: abstract, §I, §II.C, §III, §V.D, §VI.
  - Claim, abstract: all space-times under study are unstable "(the extremal case a = 2M is not analyzed separately)";
    "the leading unstable black-hole and wormhole modes approach the same limiting growth rate" M·0.10121 as a/(2M) → 1
    from either side. §II.C: "at a = 2M we obtain a regular extremal black hole with a horizon at minimum radius". §V.D:
    "This does not establish the regularity of the perturbations at the extremal horizon, which requires a separate
    analysis." §VI: "the regularity of perturbations at the degenerate horizon of the extremal configuration a = 2M
    remains to be established".
  - Gauge, §III.A: "The gauge δβ = 0 is attainable only for perturbations with δβ = 0 at the throat, whereas the gauge
    δψ = 0 is always attainable. This is the source of all subtleties discussed below."
  - Contains: the radial spectrum of the whole Simpson–Visser family except a = 2M, the rate of the Ellis limit
    (1.18197). Does not contain: a horizon at the minimum with κ ≠ 0, the horizon sector, c₂, an exact mode. Bears on S1
    (gauge), S3 (a = 2M is extremal, κ = 0; the lemma excludes it) and on S2 by contrast (the limit of the rate is not
    κ).
- **[FLV24] Franzin, Liberati & Vellucci 2024**, "From regular black holes to horizonless objects: quasi-normal modes,
  instabilities and spectroscopy", JCAP 01 (2024) 020, arXiv:2310.11990 v2, doi:10.1088/1475-7516/2024/01/020 (INSPIRE
  and arXiv checked). Read: abstract, §III, §V. Gravitational and test-field QNMs of Bardeen-like and Simpson–Visser
  black-hole mimickers with a phantom + magnetic source, up to the one-way wormhole with an "extremal null throat" at ℓ
  = 2M; no radial (monopole) sector; the Aretakis instability is mentioned as conjectured for the extremal case. Bears
  on S3 only as the extremal contrast.
- **[RJLR25] Róis, Junior, Lobo & Rodrigues 2025**, "Horizons, throats and bounces in hybrid metric-Palatini gravity
  with a non-zero potential", JCAP 07 (2025) 078, arXiv:2504.07861 v2, doi:10.1088/1475-7516/2025/07/078 (INSPIRE and
  arXiv checked). Read: abstract, §III.C, §V.C, §VI.
  - §V.C: "The former type of solutions emerges when exactly one horizon exists – either an EH or an extremal EH – and
    it is also a throat (xT = xH ), corresponding to an extremal null throat." "occurs when the minimum of R(x)
    coincides with the outermost horizon (among several), whether it is an EH or an extremal EH". §VI: "guiding future
    research into the stability properties of these solutions".
  - Contains: charged phantom-sector solutions with the minimum of R on a horizon, with Penrose diagrams. Does not
    contain: a perturbation analysis; per solution, whether the horizon is simple. Bears on the list of solutions in
    S3's scope.
- Read at the passage (G): 2409.12101 (Alencar, Nilton, Rodrigues & Silva, arXiv v2 25 Aug 2025; INSPIRE's title is
  garbled by MathML, journal not read; §VI: "the WH’s throat coincides with the event horizon, and, with this, we have a
  black throat"), 2603.05543 (Alencar, Crispim, Sáez-Chillón Gómez & Silva, Phys. Dark Univ. 53, 102427 (2026); §II.C:
  "the root of A(x) at x = 0 is degenerate (extremal horizon)"), 2609.14199 (Alencar, Bronnikov, Crispim, Sáez-Chillón
  Gómez & Silva 2026, arXiv v2 16 Sep 2026, no journal reference on INSPIRE; §II.C: "There can be an intermediate case
  where the minimum of r(x) coincides with a horizon"), 2201.01650 (Franzin, Liberati, Mazza, Dey & Chakraborty, PRD
  105, 124051 (2022); rotating Simpson–Visser, null throat on an extremal horizon, test scalar). All solution-level or
  test-field.
- **[XCY25] Xu, Chew & Yeom 2025**, "Dynamics of Bronnikov–Ellis wormhole with double-null simulation", JCAP 08 (2025)
  012, arXiv:2503.07610 v1, doi:10.1088/1475-7516/2025/08/012 (INSPIRE and arXiv checked). Read: abstract, §III.B.1,
  §IV. Nonlinear double-null evolution of the Ellis wormhole; a pulse "splits the throat into two distinct apparent
  horizons" (§III.B.1). Wormhole, no Killing horizon at the throat. Context for the nonlinear fate; not for S1–S4.

**Horizon instabilities (S3's contrast): κ = 0, test fields, inner horizons**
- **[Are11a] Aretakis 2011**, "Stability and Instability of Extreme Reissner–Nordström Black Hole Spacetimes for Linear
  Scalar Perturbations I", CMP 307, 17, arXiv:1110.2007 v1, doi:10.1007/s00220-011-1254-5; **[Are11b]** "… II", Ann.
  Henri Poincaré 12, 1491, arXiv:1110.2009 v1, doi:10.1007/s00023-011-0110-7 (INSPIRE and arXiv checked). Read: abstract
  and §1 of I; II by abstract and patterns. Abstract of I: "The fundamental new aspect of this problem is the degeneracy
  of the redshift on H+." §1.1.1: "a black hole is called extreme if the surface gravity vanishes on H+". Linear scalar
  field on extreme RN; non-decay and blow-up of transversal derivatives. Does not contain: non-extremal horizons, the
  metric perturbation.
- **[Are15] Aretakis**, "Horizon Instability of Extremal Black Holes", Adv. Theor. Math. Phys. 19, 507 (2015),
  arXiv:1206.6598 v2, doi:10.4310/ATMP.2015.v19.n3.a1 (INSPIRE and arXiv checked). Read: abstract, §1, §2.1. Abstract:
  "in stark contrast with the subextremal case for which decay is known for all derivatives along the event horizon".
  Assumption A2, Eq. (2): ∇_V V = 0 on H⁺; "The condition (2) captures the extremality of H+." The instability is
  "completely determined by local properties of extremal horizons".
- **[LMRT13] Lucietti, Murata, Reall & Tanahashi 2013**, "On the horizon instability of an extreme Reissner–Nordström
  black hole", JHEP 03 (2013) 035, arXiv:1212.2557 v2, doi:10.1007/JHEP03(2013)035 (INSPIRE and arXiv checked). Read:
  abstract, §1, contents. Abstract: the instability occurs for a massive scalar and for "coupled linearized
  gravitational and electromagnetic perturbations"; "instabilities are present for initial perturbations supported
  outside the horizon, e.g. an ingoing wavepacket". Extreme RN only. Does not contain: a non-extremal horizon, an area
  extremum.
- **[MRT13] Murata, Reall & Tanahashi 2013**, "What happens at the horizon(s) of an extreme black hole?", CQG 30,
  235007, arXiv:1307.6800 v1, doi:10.1088/0264-9381/30/23/235007 (INSPIRE and arXiv checked). Read: abstract, §1.
  Abstract: "generically the endpoint of the instability is a non-extreme RN solution". §1: "So, for a non-extreme black
  hole, outgoing radiation at the horizon decays. However, for an extreme black hole, κ = 0 so the horizon redshift
  effect is absent: outgoing radiation at the horizon does not decay."
- Read at the abstract and the patterns (G or short S): Angelopoulos, Aretakis & Gajic 1809.10037 (PRL 121, 131102
  (2018)) and 1802.05692 (CQG 35, 155007 (2018); data through the bifurcation sphere of sub-extremal RN, polynomial
  tails); Casals, Gralla & Zimmerman 1606.08505 (PRD 94, 064003 (2016), extremal Kerr); Gralla, Ravishankar & Zimmerman
  1911.11164 (JHEP 05 (2020) 094, extremal BTZ); Katagiri & Kimura 2112.09832 (PRD 105, 2022; "general static and
  spherically symmetric extremal black hole spacetimes", areal radius as the coordinate); Porfyriadis, Rosen &
  Tsaraktsidis 2512.23629 (JHEP 06 (2026) 068); Agrawal, Charalambous, Donnay, Liberati & Neri 2604.15260 (JHEP 11
  (2026) 104); Solodukhin 2610.11664 (arXiv v3, 17 Nov 2026; extremal mimickers, "Aretakis-like modes arise in a local
  analysis of the throat, they do not lead to a physical instability of the mimicker"). All have κ = 0 (or
  near-extremal) and a test field.
- **[CDLPV18] Carballo-Rubio, Di Filippo, Liberati, Pacilio & Visser 2018**, "On the viability of regular black holes",
  JHEP 07 (2018) 023, arXiv:1805.02675 v1, doi:10.1007/JHEP07(2018)023 (INSPIRE and arXiv checked). Read: §3. §3, after
  Eq. (3.6): "the mass parameter in the region A grows exponentially, on a timescale determined by the surface gravity
  of the inner horizon". An inner-horizon growth at the surface-gravity rate with generic data (ingoing tails and an
  outgoing flux), in the mass function; a different mechanism from S3.

**Gauge fixing of a horizon, minimal sections, trapping horizons (S3)**
- **[HW13] Hollands & Wald 2013**, "Stability of Black Holes and Black Branes", CMP 321, 629, arXiv:1201.0463 v4,
  doi:10.1007/s00220-012-1638-1 (INSPIRE and arXiv checked). Read: §2.1.
  - §2.1: "the linearized Raychaudhuri equation implies that, to first order in λ, the expansion must be constant along
    the horizon. Thus, the perturbed expansion δϑ must vanish on the horizon at all times"; "We now wish to argue that
    this condition can always be imposed."; with `K = κ[u∂_u − r∂_r]` (Eq. 9), "if λ1 = 0, we would have ϑ(1) = 0, i.e.
    B1 would be a MOTS lying outside the black hole, which is impossible", so "We now show that λ1 is strictly positive,
    λ1 > 0" and "we can uniquely solve the equation C(f ) = j for any smooth j". Lemma 1 is for "a stationary black hole
    spacetime satisfying the vacuum Einstein equations". "For another argument leading to the same conclusion in the
    more general context of “non-evolving horizons”, see Prop. 3 of [26]" ([Mars12]).
  - Contains: the gauge condition δϑ|_B = 0 and its proof for vacuum stationary black holes. Does not contain: matter
    with a potential, a horizon whose section is minimal, a growth rate of δr. Bears on S3 (with [Mars12]).
- **[Mars12] Mars 2012**, "Stability of MOTS in totally geodesic null horizons", CQG 29, 145019, arXiv:1205.1724 v1,
  doi:10.1088/0264-9381/29/14/145019 (INSPIRE and arXiv checked). Read: abstract, §1, Props. 2–3, Cor. 3, Thm. 1.
  - Abstract: "We link these stability properties with the surface gravity of the horizon and/or to the existence of
    minimal sections." Eq. (11): `L ψ = −Δψ + 2 s_ℓ(∇ψ) + ψ(2 div s_ℓ − ½κ_ℓθ_k + ½⟨k, tr N⟩)`; Killing horizons have N
    = 0 ("This definition includes all Killing horizons"). §1: "non-evolving horizons with non-zero and constant surface
    gravity are marginally stable if and only if they admit a section with vanishing total null expansion". Thm. 1: "If
    H admits a minimal section (i.e. a section S0 with vanishing mean curvature vector H|S0 = 0), then H is marginally
    stable."
  - Contains: the stability operator of a section of a Killing horizon and its vanishing on minimal sections, with the
    Einstein tensor in Eq. (13). Does not contain: perturbations of the spacetime, δr, a scalar source, a growth rate.
    Bears on S3.
- **[Hay09] Hayward 2009**, "Wormhole dynamics in spherical symmetry", PRD 79, 124001, arXiv:0903.5438 v1,
  doi:10.1103/PhysRevD.79.124001 (INSPIRE and arXiv checked). Read: abstract, §I, §IV. §I: "a Morris-Thorne wormhole
  throat is a double trapping horizon, which will generally bifurcate under a dynamic perturbation, such as someone
  crossing it". §IV: "In static cases, a wormhole mouth as defined above must be a double trapping horizon". Timelike
  throats; no growth rate. Also read: [Hay99] (entry above, patterns), Hayward 2002 gr-qc/0203051 and Hayward & Koyama
  2004 gr-qc/0406080 (patterns), 2104.08935, 2104.11134, 2601.05902 (patterns; no statement on S1–S4).
- Physical-process first law, read at the passage: Mishra, Chakraborty, Ghosh & Sarkar 1709.08925 (JHEP 09 (2018) 034:
  "The event horizon exhibits a teleological property and one has to impose future boundary conditions to ensure the
  stability under perturbation"), Visser & Yan 2403.07140 (JHEP 10 (2024) 029, Eq. (1.19): dδθ_v/dv = −8πG δT_vv),
  Hollands, Wald & Zhang 2402.00818 (PRD 110, 024070 (2024)). The Raychaudhuri statements agree with [HW13]; none treats
  an area extremum.

*Search limits (#260).* INSPIRE citers of [BKZ12], [BF06], [BMD07], [SV19], [Are15], [LMRT13], [MRT13] (262, 453, 197,
543, 245, 121, 105 records) downloaded whole; 15 other INSPIRE queries; 13 arXiv API queries; 8 WebSearch queries in
extended mode (`I260-mac/queries/`). 424 full texts searched with the 13 patterns of `I260-mac/screen_patterns.txt`; 59
of them read at the depth stated per entry; 52 citers of [BKZ12] without an arXiv id not read. A statement worded unlike
every pattern would be missed. AI-directed text: none (27 pattern matches, all author names, "do not report", or
disclosures of AI use).

### Other proposals of a new universe inside a black hole
- **[FMM89] Frolov, Markov & Mukhanov 1989**, "Through a black hole into a new universe?", PLB 216, 272,
  doi:10.1016/0370-2693(89)91114-3, and **[FMM90]** "Black holes as possible sources of closed and semiclosed worlds",
  PRD 41, 383, doi:10.1103/PhysRevD.41.383. No arXiv text. Bibliographic data checked on INSPIRE; the publisher pages
  answered HTTP 403, so only the INSPIRE abstracts were read.
  - Claim ([FMM90] abstract): "instead of the singularity the closed world can be formed inside the black hole". "After
    passing the deflation stage the closed world in the black-hole interior may begin to inflate and give rise to a new
    macroscopic universe." "The described model may be considered as an example of the creation of a closed or
    semiclosed world 'in the laboratory.'"
  - Mechanism (abstracts): "on the assumption that some limiting curvature exists", "the Schwarzschild metric inside the
    black hole can be attached to the de Sitter one at some spacelike junction surface which may represent a short
    transition layer", with Israel's thin-shell method. The abstracts name no matter model or energy condition. [EB01]
    says of [FMM89] and a 1992 paper of Morgan's that the matching was "postulated without any constructive realization
    (but based on the limiting curvature conjecture)".
  - Not determinable from the abstracts: formation from collapse, stability, the shell's energy conditions. Full text
    not read.
  - Reach (#224): INSPIRE gives IC/88/91 and KEK scan 2000-33-351; the KEK server rejected plain requests, the ICTP and CERN
    addresses tried failed, OSTI has the abstract only. A Schwarzschild interior glued to the expanding de Sitter wedge at
    r = r₀ needs a shell with p_x < 0 for every r₀, r_s and H (`baby_universes_junction.md` §7, *checked symbolically*, and
    rederived by the coordinator); whether [FMM90] state this is not known.
- **[Smo92] Smolin 1992**, "Did the universe evolve?", CQG 9, 173–192 (INSPIRE; [Smo06] prints 173–191),
  doi:10.1088/0264-9381/9/1/016, and **[Smo06]**, "The status of cosmological natural selection", arXiv:hep-th/0612185
  v1 (preprint, no journal reference). Read [Smo06] §1–3 in full and the INSPIRE abstract of [Smo92].
  - Claim ([Smo92] abstract): "all final singularities 'bounce' or tunnel to initial singularities of new universes at
    which point the dimensionless parameters of the standard models of particle physics and cosmology undergo small
    random changes". [Smo06], hypothesis I: "Black hole singularities bounce and evolve to initial states of expanding
    universes". The "master prediction": "almost no change in the parameters of the standard model from the present
    values will increase the numbers of black holes produced."
  - Mechanism: none stated. The bounce is cited to Wheeler and, as "theoretical evidence", to Bojowald's 2005 Living
    Reviews article (ref. 14; [Smo06] does not print its title). No matter, no energy condition, no solution. The
    hypothesis is about an ensemble of universes and its falsifiable consequences, listed in [Smo06] §3 (a neutron-star
    mass limit, single-parameter inflation, little early star formation).
  - Does not contain: any geometry of the bounce, its formation or its stability. [Smo06] §4, a reply to Vilenkin, was
    not read.
- **[EB01] Easson & Brandenberger 2001**, "Universe generation from black hole interiors", JHEP 06 (2001) 024,
  arXiv:hep-th/0103019 v1, doi:10.1088/1126-6708/2001/06/024 (INSPIRE and arXiv checked). Read: abstract, §I, start of
  §II.
  - Claim: "scenarios in which the universe is born from the interior of a black hole may not posses [sic] many of the
    problems of the Standard Big-Bang (SBB) model", the horizon, flatness and structure problems, "not necessarily
    requiring a long period of cosmological inflation". "Our conclusions are completely independent of the details of
    general models."
  - Mechanism: a Schwarzschild–de Sitter matching (Fig. 1), motivated by limiting curvature; the paper reviews the
    realizations it knows. The argument for horizon and homogeneity is that the forward light cone of a horizon point
    covers the radial direction of the de Sitter universe and that a particle may travel "about 180°" round the black
    hole. No energy condition is discussed (searched). Does not contain: a specific solution, formation, stability.
- **[BHR21] Brandenberger, Heisenberg & Robnik 2021**, "Through a black hole into a new universe", IJMPD 30, 2142001,
  arXiv:2105.07166 v1, doi:10.1142/S0218271821420013 (INSPIRE and arXiv checked). Read in full (an essay).
  - Claim: "an S-Brane which arises in the inside of the black hole horizon when the Weyl curvature reaches the string
    scale induces a continuous transition between the inside of the black hole and the beginning of a new universe."
  - Mechanism: a spacelike brane with "positive tension but vanishing energy density. This object violates the usual
    energy conditions". Israel matching of the Schwarzschild interior to a homogeneous anisotropic cosmology. "Making
    use of the equations of General Relativity coupled to matter obeying the usual energy conditions is not possible".
    Late-time isotropy needs "matter which acts for some time as an effective cosmological constant", which "is very
    difficult to obtain" from string theory in effective field theory.
  - Does not contain: a collapse that forms the brane ("We assume that there is no matter present to the past of the
    S-brane"), stability, details (deferred to arXiv:2103.02842, not read).
  - Junction check (#224, `baby_universes_junction.md`): (4)–(5) are the Israel conditions of the Schwarzschild interior and a Kantowski–Sachs cosmology with √g read as √(−f(r₀)); the printed right-hand side S is twice the one that the action (1) gives at κ = 1. Positive tension is needed for ḃ/b > 0 at every r₀ < r_s and for the volume at r₀ < 3r_s/4; the brane violates the NEC and WEC for σ > 0.
- **[Pop10] Popławski 2010**, "Radial motion into an Einstein–Rosen bridge", PLB 687, 110, arXiv:0902.1994 v3,
  doi:10.1016/j.physletb.2010.03.029, and **[Pop16] Popławski 2016**, "Universe in a black hole in Einstein–Cartan
  gravity", ApJ 832, 96, arXiv:1410.3881 v2, doi:10.3847/0004-637X/832/2/96 (INSPIRE and arXiv checked). Read:
  abstracts, introductions, closing paragraphs.
  - Claim ([Pop10]): "observed astrophysical black holes may be Einstein-Rosen bridges, each with a new universe inside
    that formed simultaneously with the black hole. Accordingly, our own Universe may be the interior of a black hole
    existing inside another universe." [Pop16]: "The collapsing matter in a black hole should therefore bounce at a
    finite density and then expand into a new region of space on the other side of the event horizon, which may be
    regarded as a nonsingular, closed universe."
  - Mechanism: [Pop16] spin–torsion coupling "generates gravitational repulsion at extremely high densities in fermionic
    matter, approximated as a spin fluid", plus particle production for a finite inflation. It does not name an energy
    condition. [Pop10] notes that the Einstein–Rosen bridge metric needs a source that "violates the energy conditions
    at the throat".
  - Does not contain: a collapse computation ([Pop16]: "The formation of trapped surfaces and wormhole throats ... could
    be explored numerically in the ECSK theory"; the merging of throats is conjectured), stability of the bounce.
- **[Hay99] Hayward 1999**, "Dynamic wormholes", IJMPD 8, 373, arXiv:gr-qc/9805019 v3, doi:10.1142/S0218271899000286
  (INSPIRE and arXiv checked). Read: abstract, first section. The framework [SH02] tests: "the null energy condition is
  violated everywhere on a generic wormhole horizon"; "if a wormhole's negative-energy source fails, it may collapse
  into a black hole. Conversely, irradiating a black-hole horizon with negative energy could convert it into a wormhole
  horizon." No universe claim.
- **[Hay06] Hayward 2006**, "Formation and evaporation of regular black holes" (arXiv title "non-singular"), PRL 96,
  031103, arXiv:gr-qc/0506126 v2, doi:10.1103/PhysRevLett.96.031103 (INSPIRE and arXiv checked). Read: abstract,
  introduction. Gives regular space-times that "describe the formation of a (locally defined) black hole from an initial
  vacuum region, its quiescence as a static region, and its subsequent evaporation". The static region "behaving as a
  cosmological constant at small radius", with a regular centre; infalling positive-energy flux in collapse and
  negative-energy flux in evaporation. A formation model of a regular black hole, with no new-universe claim and no
  throat.
- **[MSOF23] Masó-Ferrando, Sanchis-Gual, Font & Olmo 2023**, "Birth of baby universes from gravitational collapse in a
  modified-gravity scenario", JCAP 06 (2023) 028, arXiv:2304.12018 v2, doi:10.1088/1475-7516/2023/06/028 (INSPIRE and
  arXiv checked). Read: abstract, §4 results, §5 conclusions.
  - Claim: in the Palatini f(R) = R + ξR² frame, collapsing boson stars show "the formation of a finite-size,
    exponentially-expanding baby universe connected with the outer (parent) universe via a minimal area surface (a
    throat or umbilical cord)", "at all times hidden inside a horizon". In the Einstein frame "the endpoint is a
    nonrotating black hole surrounded by a quasi-stationary cloud".
  - Mechanism: a complex canonical scalar in a metric-affine theory, mapped to an Einstein-frame problem; the baby
    universe is read in the f(R)-frame metric. The paper names no energy condition in the passages read.
  - Contains: a numerical relativity evolution (BSSN, unstable boson-star branch): horizon at t ≈ 80, throat at t ≈
    84.6, apparent horizon at t ≈ 87; the expansion is "exponential and superluminal" for t ∈ [84.6, 91.8]. Does not
    contain: a run past t = 91.8 ("the region between the center and the maximal sphere becomes unreachable", from the
    slicing); a stable throat (it expects the throat to close as the scalar cloud is absorbed); the baby universe's
    matter content.
- **[BCHM25] Bueno, Cano, Hennigar & Murcia 2025**, "Dynamical formation of regular black holes", PRL 134, 181401,
  arXiv:2412.02742 v2, doi:10.1103/PhysRevLett.134.181401 (INSPIRE and arXiv checked). Read: abstract and passages on
  the new universe.
  - Claim: in a theory with "an infinite tower of higher-derivative corrections", "the shell bounces into a new universe
    through a white hole explosion"; the thin dust shell collapse is solved analytically and is "singularity-free,
    geodesically complete".
  - Mechanism: higher-curvature gravity with second-order equations on spherical metrics, whose terms give the bounce;
    the matter is a thin dust shell. No energy condition was read. Contains: "stable time evolution within spherical
    symmetry" of the theory, the shell collapse, horizons above a critical mass. Does not contain: a scalar or phantom
    source, or the stability of the inner horizon ("can also be studied in the spherically symmetric sector", left for
    later).

### Named, not read
- Farhi, Guth & Guven, "Is it possible to create a universe in the laboratory by quantum tunneling?", NPB 339, 417
  (1990), doi:10.1016/0550-3213(90)90357-J. Oshita & Yokoyama, "Creation of an inflationary universe out of a black
  hole", PLB 785, 197 (2018), arXiv:1601.03929. Chakrabarty, Abdujabbarov, Malafarina & Bambi, "A toy model for a baby
  universe inside a black hole", EPJC 80, 373 (2020), arXiv:1909.07129. Dymnikova, "Universes inside a black hole with
  the de Sitter interior", Universe 5, 111 (2019). Named in [BHR21]'s reference [3]. The INSPIRE records match; the
  texts and abstracts were not read.
- Titles only (INSPIRE citers of [SV19]): Lobo et al., "Dynamic thin-shell black-bounce traversable wormholes", PRD
  2020, arXiv:2003.09419; Franzin et al., arXiv:2201.01650, 2207.08864 and 2310.11990 (perturbations and stability of
  regular black holes and wormholes, rotating and not); Bueno et al. as above. [FMM89] and [FMM90] full texts. [BHR21]'s
  companion arXiv:2103.02842.

*Search limits.* The papers beyond the issue's list came from [BHR21]'s reference [3] and from INSPIRE citer lists of
[BF06] and [SV19] (`refersto:recid:` queries on stability, collapse, formation and dynamics in the title; raw
`inspire_citer_searches.txt`). No other search was made. Papers that cite [SH02] were covered by #50
(`critical_ghost.md` §6). AI-directed text: none found.

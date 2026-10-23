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

## Cosmology / observation
- Farnes 2018, A&A 620, A92; critique Socas-Navarro 2019, A&A 626, A5 (large-scale runaway).
- Trivedi & Loeb 2026, arXiv:2605.10976 — GW signatures of negative-mass binaries.
- ALPHA-g 2023: antihydrogen falls down, g_H̄ ≈ 0.75 ± 0.29 g; ALPHA-g, AEgIS, GBAR target 1%.

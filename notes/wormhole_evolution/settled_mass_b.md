# The settled mass of family G against the normal amplitude B, and the static member at departure (#266)

Code: [settled_mass_b.py](settled_mass_b.py), a thin driver that imports #41's `settled_mass.py`, #253's `settled_mass_g.py`, #40's `threshold.py`, #225's `throat_locus.py`
and #39's `dn.py` unchanged. Output of record: [settled_mass_b_output.txt](settled_mass_b_output.txt) (`python3 settled_mass_b.py analyze --raw=<negative-mass-raw>/I266-mac`, exit 1:
one gating check fails, REL, below). Raw output, with the criteria written before the runs of record: negative-mass-raw `I266-mac` (README with the criteria, amendment A1 and the
mutation targets; `runs/`, `mut/`, `dev/`, `scripts/`). Status words as in AGENTS.md: every mass and every member parameter below is *numerical evidence*; what is not shown is
marked *conjecture*. No sentence here says a result is new. Units: a = 1 is the throat radius of the initial static wormhole, 8πG = 1, D = 4, M = r/2 on the apparent horizon.

## Answer

**B = 0 reproduces S's settled mass, B up to 0.1 gives a settled mass that rises with the asymmetry of the member the run leaves, and B = 0.2 has no stall to measure.** At D = 4,
δ = 1e-8, h = 0.005, the plateau of M_AH to e^{−f} = 1.57e6 on the horizon is 0.223086 a at B = 0 (S: 0.22309, difference −3.8e-6; bar 6e-4), 0.223563 a at B = 0.025, 0.225235 a at
B = 0.05 and 0.235548 a at B = 0.1 (#253's 0.235548). The member the run leaves, read where the unstable mode has grown to 1e-3 of the throat radius, is r_th = 1.0000, 0.99945, 0.99836,
1.00288 a with γ₁ = 0, 0.00294, 0.01184, 0.04832: γ₁ ≈ 4.8 B² (local exponent 2.0). The asymmetry rises as B², the radius first dips and then rises, and the mass follows.

At B = 0.2 the near-threshold run does not stall: T = 5.3 on the collapse side and 12.2 on the expansion side for every A within 1e-3 of p*ₕ, and T does not grow as p*ₕ is approached.
The chain stops at "horizon lost" (e^{−f} = 0.77 on the horizon, 113 columns at h = 0.005); the first trapped cell has r = 1.112 a and M = 0.556 a. There is no member at departure and
no settled mass at B = 0.2: this is a finding about G, not a failed check. The stall ends at B ≈ 0.16 (section "Where the stall ends").

| acceptance (issue #266) | verdict | numbers |
|---|---|---|
| B = 0 reproduces 0.22309 a within 6e-4 | **met** (RB0) | G(B = 0) and O, bit for bit the same: 0.223086 (h = 0.005), 0.223095 (0.01), 0.223143 (0.02); differences −3.8e-6, +5.3e-6, +5.3e-5; Richardson 0.223083 |
| M(B) at two h, against γ₁ and r_th at departure | M(B) at two h **met**; against the member: table below, one gating check **not met** | |plat(0.01) − plat(0.005)| ≤ 1.5e-5 for every B (RH, bar 5e-5); γ₁, r_th at the two h agree to 1e-6 (HST); REL below |
| second order | **met, post hoc** (RC6-ph, PO-ph) | needs three h; h = 0.02 is post hoc by the coordinator's ruling. M_AH(x⁺) orders 1.71 to 2.34 in 18 windows (band 1.5 to 2.5); p*ₕ order 2.000 for every stalling member |
| a planted member through the same code, or the reason it cannot be | **met** (PMEM) for the extraction; the chain **cannot** be started from one | exact static members, γ₁ = 0 to 0.3, recovered to 1e-6 (γ₁), 1.4e-8 (r_th), 1.6e-6 (m±); reason below |
| each check has a `--mutate=` whose effect is computed before it is declared | **met** | 19 mutations, targets computed on the development set (raw README); on the runs of record: see "Controls" |

Controls: 37 checks, 20 gating. **One gating check fails: REL** (7.9e-5 against 6e-5, below). The 17 post hoc checks all pass. The criteria were not moved; amendment A1 of the raw README (the split
into gating and post hoc checks, the coordinator's ruling) was committed before any run of record was analysed.

## What was run

- Eight members: O, G with B = 0, 0.025, 0.05, 0.1, 0.2 (the issue's), and G with B = 0.125 (post hoc, added by this note), at h = 0.01 and 0.005 (the issue's) and 0.02 (post hoc). 21 jobs,
  each: p*ₕ bisected to hi − lo ≤ 1e-14 (#40's `threshold.bisect`, 43 to 45 runs, 5 s at h = 0.02 to 75 s at 0.005), #41's chain at p*ₕ − 1e-8 (1 to 6 s, to e^{−f} = 1.54e6 to 1.57e6), and the
  pair p*ₕ ∓ 1e-8 through `throat_locus.sweep` (7 to 20 s). All exit codes 0. The whole set took 15 min on two lanes (nice 10, one thread each). "1e-14 relative" in the issue is read as #40's
  hi − lo ≤ 1e-14, reached to 5.7e-15 to 9.1e-15 (1e-9 of p*ₕ at B = 0, 1e-13 at B = 0.2, 1e-6 of δ).
- **My p*ₕ at B = 0.1 equals #40's** (dh610) to 4.5e-14, 1.7e-14, 4.0e-14 at h = 0.02, 0.01, 0.005 (bar 1e-10, the h of its job.log equal to the h passed). G with B = 0 and O have the same p*ₕ bit for bit
  at every h, and the same chain series: a zero normal pulse changes nothing.
- The chain is #253's driver unchanged; B = 0.1 reproduces #253's plateau means at h = 0.02, 0.01, 0.005 to 1.5e-5, 1.1e-6, 7.5e-8 (RB10, bar 1e-4). The first law with both fields holds to 5.2e-4 (RC4).
- **The member at departure** (#225's method, `slice_masses` and `fit_asym` on profile slices every 0.5 in t): the pair p*ₕ ∓ 1e-8 is the run of record's amplitude and its mirror. Its mean profile has
  no first-order unstable mode, so the member is read on it: the throat radius is the mean of the two locus minima, s_th is #225's fit on the mean profile (radii 1.25, 1.5, 2, 3, 4 r_th), γ₁ = s_th/√(1 − s_th²).
  Departure is the first t at which the mode's part |r_min(p*ₕ + δ) − r_min(p*ₕ − δ)|/2 reaches ε = 1e-3 of r_th (t = 13.29 at B = 0, 13.30, 13.36, 13.61, 13.86 at B = 0.025, 0.05, 0.1, 0.125).
  ε and the choice of profile are mine (the issue gives no definition); both are tested below.

## M(B) against the member at departure

h = 0.005, δ = 1e-8; M the plateau mean over the last 1.0 of x⁺; dM′ = M/r_th − M(B = 0). The B = 0.125 row is post hoc.

| B | p*ₕ (h = 0.005) | γ₁ | r_th | M | M − M₀ | M/r_th | dM′ | dM′/γ₁ | m₊ of the member |
|---|---|---|---|---|---|---|---|---|---|
| O, 0 | 5.9423e-6 | 0 (|γ₁| ≤ 3e-8) | 1.00000 | 0.223086 | 0 | 0.223087 | 0 | | 0 |
| 0.025 | 8.8746e-4 | 0.00294 | 0.99945 | 0.223563 | +0.00048 | 0.223687 | +0.00060 | 0.205 | +0.00292 |
| 0.05 | 3.5695e-3 | 0.01184 | 0.99836 | 0.225235 | +0.00215 | 0.225606 | +0.00252 | 0.213 | +0.01161 |
| 0.1 | 1.4860e-2 | 0.04832 | 1.00288 | 0.235548 | +0.01246 | 0.234872 | +0.01179 | 0.244 | +0.04497 |
| 0.125 (post hoc) | 2.3906e-2 | 0.07545 | 1.01588 | 0.247100 | +0.02401 | 0.243237 | +0.02015 | 0.267 | +0.06828 |
| 0.2 | 6.8270e-2 | no stall | | M at the first trapped cell 0.556 | | | | | |

- **The asymmetry is second order in B**: γ₁/B² = 4.70, 4.74, 4.83, 4.83 (0.025 to 0.125); the tuned ghost amplitude is too: p*ₕ(h→0)/B² = 1.41, 1.43, 1.49, 1.53 (Richardson, three h; p*∞ = 4e-10 at B = 0).
- **The radius is not monotone**: r_th − 1 = −5.5e-4, −1.6e-3, +2.9e-3, +1.6e-2: it dips to −1.6e-3 at B = 0.05 and then rises. The mass is not a function of r_th alone.
- The local exponents d ln/d ln B between consecutive B: γ₁ 2.0; dM′ 2.07, 2.23, 2.40; M − M₀ 2.17, 2.54, 2.94. The mass shift grows faster than γ₁ because dM′/γ₁ is not constant.
- **A relation along B** (REL): dM′ = a γ₁ + b γ₁², fitted to the other points of the set and predicted for the one left out ((a, b) = (0.2027, 0.853) by least squares through the four; the leave-one-out fits give 0.2026 to 0.2027 and 0.853 to 0.856):

| held out | measured dM′ | predicted | error | the same without the r_th scaling (M − M₀) |
|---|---|---|---|---|
| 0.025 | +0.00060 | +0.00060 | +1.4e-6 | −1.2e-5 |
| 0.05 | +0.00252 | +0.00252 | +6.7e-7 | −5.3e-5 |
| 0.1 | +0.01179 | +0.01178 | −2.0e-6 | +9.6e-5 |
| 0.125 | +0.02015 | +0.02016 | +5.5e-6 | −2.6e-4 |

  Over four B (post hoc, `REL-ph`) the largest error is 5.5e-6 against the bar 6e-5. **The gating REL, over B = 0.025, 0.05, 0.1 only, is not met**: 7.9e-5 against 6e-5 (held out B = 0.1,
  predicted from the two smaller B, 4 times smaller in γ₁); 1.2e-6 and 3.8e-6 at the other two. The bar was declared for the four-point version and was not moved. The scaling by r_th matters:
  without it the four-point errors are 10 to 50 times larger (up to 2.6e-4).
- **What this does not show.** The form was chosen after the development runs at h = 0.02 and 0.01 on these members, so the record at h = 0.005 reproduces it and does not test it blind. It is a
  relation along a one-parameter family: B sets γ₁, r_th, the tuned ghost amplitude and the remnants of the pulses together, and any smooth M(B) is a smooth function of γ₁(B). So M following
  the member, and not something else that also scales with B², is still a *conjecture*. #267 varies δ at fixed B, which moves the member along the stall and not along B.
- The member's end masses (m₊, m₋ = +0.0450, −0.0523 at B = 0.1) are #225's 0.045, −0.052. M − M₀ is 0.16, 0.19, 0.28 of m₊ at B = 0.025, 0.05, 0.1: not a fixed fraction.

### Sensitivity of the member to its definition

| B | h | ε = 1e-3: t_dep, r_th, γ₁ | ε = 1e-2: t_dep, r_th, γ₁ | collapse side alone, ε = 1e-3: r_th, γ₁ |
|---|---|---|---|---|
| 0 | 0.005 | 13.287, 1.00000, 0 | 15.234, 0.99986, 0 | 0.99896, 0 |
| 0.05 | 0.005 | 13.355, 0.99836, 0.01184 | 15.311, 0.99775, 0.01184 | 0.99733, 0.01185 |
| 0.1 | 0.005 | 13.615, 1.00288, 0.04832 | 15.618, 1.00334, 0.04816 | 1.00184, 0.04835 |
| 0.125 | 0.005 | 13.864, 1.01588, 0.07545 | 15.929, 1.01966, 0.07488 | 1.01484, 0.07550 |

From ε = 1e-3 to 1e-2 γ₁ moves by 3e-4 relative (B = 0.1, 1.6e-4 absolute) and r_th by 4.6e-4. **The relation REL depends on that choice at the 1e-4 level**: with the member read at ε = 1e-2 the
four-point leave-one-out errors are −2.6e-5, −2.0e-5, +4.9e-5, −1.4e-4 (bar 6e-5) and the coefficients are (a, b) = (0.217, 0.54) instead of (0.203, 0.85). The member drifts slowly during the
departure (r_th by 4e-4 over the two e-folds of the mode), and the relation to 6e-6 is a statement about the reading at ε = 1e-3, not a property of the member independent of when it is read. The collapse side alone
(#225's reading, here at 1e-8 where the mode has grown) has r_th lower by 1e-3 at B = 0 already, which is why the mean profile is used. At B = 0.1 the member reproduces #225's table (stall pair at
1e-9) at t = 12, 13, 14: r_th 1.00316, 1.00288, 1.00292 (#225: 1.0031, 1.0028, 1.0028), s_th 0.04841, 0.04831, 0.04823 (0.0484, 0.0483, 0.0482), to 1.2e-4 and 3e-5 (C225).

## Error bars of M (kept separate, not mixed)

| kind | what | size | what it is |
|---|---|---|---|
| value | plateau mean over the last 1.0 of x⁺, h = 0.005, δ = 1e-8 | table above | |
| h | \|plat(0.01) − plat(0.005)\| | 9.1e-6 (B = 0), 1.2e-5, 1.5e-5, 9.3e-6 (0.1), 6.7e-6 (0.125) | difference between two resolutions of this code; at second order the Richardson correction to the finest is 1/3 of it |
| roundoff realization | #253's late-window scatter of M_AH between ulp-noise realizations | 6.6e-6 at h = 0.005, 2.5e-5 at 0.01 | **quoted from #253, not remeasured here**; the same kind of number as the machine-to-machine 2e-5 of #41 |
| δ | physical, not an error | +5.7e-4 from δ = 1e-6 to 1e-8 at B = 0.1 (#253) | not measured at the other B |

The M(B) entries are quoted to 6 digits and are good to the sum of the first two, 1.3e-5 to 2.1e-5. The relation REL is tested at 6e-5.

## Where the stall ends (post hoc, gates nothing)

Lifetime T at p*ₕ ∓ 10^−e for e = 6 and 8, each p*ₕ bisected to 1e-12; a stall gains 4.4 (collapse side) and 5.5 (expansion side) from e = 6 to 8, prompt gains 0.00 to 0.04. `stallscan`.

| h | last B with a stall | first B without |
|---|---|---|
| 0.02 | 0.155 | 0.16 |
| 0.01 | 0.16 | 0.165 |
| 0.005 | 0.155 | 0.16 |

The scan runs B = 0.15 to 0.175 in steps of 0.005 at h = 0.02 and 0.01, and at h = 0.005 over 0.15, 0.155 (stall) and 0.16 to 0.175 (prompt). **The runs bracket the end of the stall between B = 0.155 and
0.165, and not more tightly**: the three h put it in (0.155, 0.16), (0.16, 0.165) and (0.155, 0.16), so it is not the same at every h and no trend in h is read from a 0.005 step. It is not extrapolated to h → 0.
The development scan at h = 0.02 over B = 0.1, 0.125, 0.15, 0.175, 0.2 had the same pattern. At B = 0.15 the bisection of the pipeline to 1e-14 ends at 9e-14 to 2e-13 ("static through X", X = 28; development)
and the chain of #41 with dxp = 12 ends by "columns exhausted" at e^{−f} ≈ 9e5, so the stall regime reached by the chain to its limit is B ≤ 0.125. Above the end the threshold is prompt: the collapse side forms
a horizon at T = 5.3 with M = 0.556 a, r = 1.11 a (B = 0.2), not a horizon that grows out of the throat after a stall, and the expansion side is at T = 12 to 14.7 with no growth toward p*ₕ. That is a different
kind of threshold from the type I stall (finite lifetime on both sides at p*ₕ ± 1e-8), and is not studied here.

## The planted member, and why the chain cannot start from one

- **Through the extraction (PMEM).** Exact static ghost wormholes of GGS09a's family with γ₁ = 0, 0.00293, 0.0483, 0.0755, 0.3 (r_th, s_th and end masses m± from the closed forms) are laid out as
  the node values of a profile slice of the sweep (r, f and the centred differences r₊, r₋ the sweep forms) at t = 13 and h = 0.01, 0.005, and go through `slice_masses`, `fit_asym`, the parabolic throat
  minimum and `asym_masses`. Errors: γ₁ ≤ 9.7e-7, r_th ≤ 1.4e-8, m± ≤ 1.6e-6 (bars 1e-5, 1e-6, 1e-4). A first version with a coarse table of r(l) gave 3e-4 to 2e-3; the table was refined before the criteria were written.
- **Through the chain: it cannot be started.** A static member has no trapped surface, so there is no first trapped cell and no horizon for the chain to follow; the horizon comes from the unstable mode,
  and the settled mass of that horizon is the unknown. On the grid the static wormhole as data (A = 0, B = 0) is not even static: it collapses at T = 7.6, 8.8, 10.0 for h = 0.02, 0.01, 0.005 (truncation
  seeds the mode) and `run_real_g` refuses it (fates −1, −1 at p*ₕ ∓ δ). The chain's planted truths stay #41's PL1 (exact Schwarzschild) and PL2 (exact ghost solution), run through the same code in #253;
  RB0 and RB10 repeat its reproduction of S's 0.22309 and of #253's 0.2355. A member as characteristic data (the initial rays of the asymmetric static wormhole) is #254's, not built here.
- Planted by symmetry: S and O are mirror symmetric, *proved* in #213, so γ₁ = 0 and r_th = 1 at B = 0. Measured: |γ₁| ≤ 2.2e-7, |r_th − 1| ≤ 3.1e-6 at the three h (MB0, bars 2e-5, 1e-4).

## Controls and mutations

All 20 gating checks pass except REL; the 17 post hoc checks pass. Each of the 19 mutations is run on the runs of record (`mut/` in the raw folder) and its fired set compared with the one declared in the raw
README (computed on the development set): **all 19 fire exactly their declared set** (MATCH). `eps_h` ran twice: its first run on the records stopped in a report line added to the analysis while the mutation runs were going (a NaN member at ε × 100 in the
ε = 1e-2 relation print); the line was guarded, `mut/eps_h.log` (crashed) and `mut/eps_h.2.log` (MATCH) are both kept. Every set contains REL when REL fails without a mutation, as the raw README says.

## What this says, and does not

- *Shown (numerical evidence):* at D = 4, δ = 1e-8, the settled mass of G equals S's at B = 0 (to 4e-6) and rises with B as M − M₀ = 0.00048, 0.0022, 0.0125 a at B = 0.025, 0.05, 0.1 (0.0240 a at 0.125);
  the member at departure is the symmetric one at B = 0 and has γ₁ = 4.8 B² for B ≤ 0.125; along B the combination M/r_th − M₀ is the quadratic 0.2027 γ₁ + 0.853 γ₁² to 6e-6 (four B, leave one out) with the member read at ε = 1e-3, and to 1.4e-4 with it read at ε = 1e-2.
- *Not shown:* that M depends on the member alone (B also sets the tuned ghost amplitude and what the pulses leave); anything at B > 0.125 in the stall regime (the chain with dxp = 12 does not reach its limit
  at B = 0.15); the settled mass at B = 0.2, where the near-threshold run is prompt and #41's chain stops at the formation of the horizon; the dependence on δ at B ≠ 0.1; D = 5, 6 (#246).
- The gating REL is not met on three points; the four-point version passes. Both are reported.

## Stale or conflicting text found (for the coordinator; none was edited)

1. `settled_mass_g.md`, recommended task 1 and `throat_locus.md` task 1: B = 0.2 is listed. At D = 4 it has no stall and the chain stops; the task's acceptance ("M(B) … against γ₁ and r_th at departure") cannot be met there.
2. `throat_locus.md`, member table: the column "s_th (γ₁)" prints s_th (0.0485, 0.0484, 0.0483, 0.0482 at t = 11 to 14), while its text gives "γ₁ = 0.0485, 0.0484, 0.0483 at t = 12, 13, 14". γ₁ = s_th/√(1 − s_th²) is 0.12 % larger:
   0.04857, 0.04847, 0.04838, 0.04834 on the collapse side at h = 0.01. Immaterial at 3 digits; this note compares s_th with the table.
3. Paper V (`paper/5-ellis-threshold/main.tex`), the threshold bullet "Whether it is a function of the member at departure is open (#266, #267)": the part of it that varies B is now measured (below).
4. `critical_ghost.md` §1.3 "Update from #253" ends "That it follows the member's parameters (γ₁, r_th at departure) is a *conjecture*": still true; the new relation is along B only.

## Proposed text (for the coordinator; none applied)

**`settled_mass_g.md`**, "What this says, and does not", the conjecture bullet: append "Tested along B in #266 (`settled_mass_b.md`): M/r_th − M(B = 0) = 0.2027 γ₁ + 0.853 γ₁² for B = 0.025 to 0.125 to 6e-6, with γ₁ = 4.8 B² the
asymmetry of the member at departure; a relation along one family, not a test that M depends on the member alone." Recommended task 1: "done in #266; B = 0.2 has no stall."

**`critical_ghost.md` §1.3**, after the "Update from #253" bullet:

> *Update from #266* (`wormhole_evolution/settled_mass_b.md`; *numerical evidence*). With the normal amplitude B of family G varied at D = 4 and δ = 1e-8, the settled mass is 0.2231 a at B = 0 (S's value, to 4e-6),
> 0.2236, 0.2252, 0.2355 a at B = 0.025, 0.05, 0.1 and 0.2471 a at 0.125. The member the run leaves has γ₁ = 0, 0.0029, 0.0118, 0.0483, 0.0755 (γ₁ ≈ 4.8 B²) and r_th = 1, 0.9995, 0.9984, 1.0029, 1.0159 a, and
> M/r_th − 0.2231 a = 0.2027 γ₁ + 0.853 γ₁² to 6e-6 over the four B (the member read where its unstable mode has grown to 1e-3 of r_th; to 1.4e-4 at 1e-2). That is a relation along one family; that the mass depends on the member alone is still a *conjecture*. The stall itself ends between B = 0.155 and 0.165
> (not the same at h = 0.02, 0.01, 0.005); above that the threshold is prompt (collapse at T = 5.3 with M = 0.556 a, expansion at T = 12), and there is no member to read.

**Paper V**, the `\pending{\#266, \#267}` paragraph at line 456: replace by a sentence and keep #267 pending:

> With the normal amplitude $B$ of family G varied, the settled mass at $\delta=10^{-8}$ is $0.2231a$ at $B=0$, the mirrored family's value, and $0.2236a$, $0.2252a$, $0.2355a$ at $B=0.025$, $0.05$, $0.1$. The member at departure
> has $\gamma_1=0.0029$, $0.0118$, $0.0483$ ($\gamma_1\approx4.8B^2$) and $r_{\rm th}=0.9995a$, $0.9984a$, $1.0029a$, and $M/r_{\rm th}-0.2231a=0.2027\gamma_1+0.853\gamma_1^2$ to $6\times10^{-6}$ over these three and $B=0.125$, with the member read where its unstable mode has grown to $10^{-3}$ of $r_{\rm th}$ (`notes/wormhole_evolution/settled_mass_b.md`).
> This is a relation along one family. Above $B\approx0.16$ the threshold has no stall and the chain does not apply.

`\pending{\#267}{…as δ→0…}` stays. `paper_numbers.py` would check 0.223086, 0.235548, 0.04832, 1.00288 and 0.2027 from `settled_mass_b_output.txt`.

## Limits

- D = 4; δ = 1e-8; h = 0.01 and 0.005 (and 0.02 post hoc); the member read at ε = 1e-3, definitions mine.
- p*ₕ is bisected on this machine (#40's at B = 0.1 reproduced to 5e-14); the roundoff-realization scatter of M_AH is #253's, not remeasured; the stall scan is at three h and its end is not extrapolated.
- B = 0.125 is added by this note and gates nothing; B = 0.15 is development only (the chain's dxp = 12 and X = 28 are #41's and #40's declared values).
- The first-order unstable mode is removed from the member by averaging the pair; the mean's second-order contamination is (1e-3)² and was not separated.
- The planted member goes through the extraction, not through the evolution or the chain.

## Recommended tasks

1. *A second curve in (γ₁, r_th).* What: G with the normal pulse's width or centre changed at tuned amplitude so that γ₁ ≈ 0.012, 0.048, 0.075 are reached by another route, and the same chain and member at h = 0.01, 0.005. Why: B moves γ₁,
   r_th and the tuned ghost amplitude together; a second route separates them. Acceptance (can fail): M/r_th − M₀ lies on 0.2027 γ₁ + 0.853 γ₁² within 6e-5, or the failure is reported. Depends on: this note, #267.
2. *The chain beyond B = 0.125.* What: B = 0.15, 0.155 with dxp enlarged (so the chain reaches e^{−f} = 1.57e6) and X = 40, h = 0.01, 0.005. Why: the quadratic is extrapolated over a factor 1.4 in γ₁ there. Acceptance: REL's prediction
   for B = 0.15 (γ₁ ≈ 0.107) within 6e-5. Depends on: this note.
3. *The end of the stall.* What: bisect B at which T(p*ₕ ∓ 1e-8) − T(p*ₕ ∓ 1e-6) falls below 2.0, at h = 0.02, 0.01, 0.005, 0.0025. Why: the scan puts the end in (0.155, 0.16), (0.16, 0.165) and (0.155, 0.16) at the three h.
   Acceptance: B_c(h) bracketed to 1e-3 at four h and its trend stated. Depends on: this note.
4. *The mass of the prompt branch.* What: a horizon finder for e^{−f} < 1 on the horizon (B ≥ 0.165), then the settled mass of the prompt collapse at B = 0.175, 0.2. Why: #41's chain stops at "horizon lost" there. Acceptance: M_AH(x⁺)
   second order at three h and a planted Schwarzschild in the same gauge. Depends on: this note; new method code.

## Files

`settled_mass_b.py` (commands `pstar`, `job`, `member`, `stallscan`, `planted`, `analyze`; `--mutate=`), `settled_mass_b_output.txt`, this note; raw `negative-mass-raw/I266-mac` (`README.md` with the criteria, amendment A1 and
the mutation targets; `runs/<member>/h<h>/{pstar,chain,member,logs}`, `runs/stallscan*`, `mut/` the 19 mutated analyses of record, `dev/` development runs and their logs, `scripts/`).

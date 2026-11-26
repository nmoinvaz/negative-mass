# Residual stage of Theorem E: the explicit normal and phantom scalars

November 2026, issue [#171](https://github.com/nmoinvaz/negative-mass/issues/171). Script
[`planar_residual.py`](planar_residual.py), its stdout [`planar_residual_output.txt`](planar_residual_output.txt). It
checks the sentence "*by hand; the residual stage with these fields is not computed here*" of
[`planar_closure.md`](planar_closure.md) §3 (#161). Conventions are those of [`two_scalars.md`](two_scalars.md)
(Setup); the planted metrics are #78's, from [`same_sign_scalars.md`](same_sign_scalars.md). The script imports
nothing from the repository.

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically** (sympy, exact, with controls that fail).

## Answer

**With φ₁ = α(cH − (4/c)J) normal, φ₂ = α(cH + (4/c)J) phantom, α = 1/(4√(2π)) = 1/(2√(8π)) and V = 0, G_ab − 8π(T₁ −
T₂)_ab vanishes exactly for the metric g = −2e^{2H(p)J(q)} dp dq + dy² + dz² with H, J undefined functions and c a
symbol, so for every c ≠ 0** (**checked symbolically**). The normalizations and V = 0 are not a choice: the residual
is linear in (α₁², α₂², V), and its unique solution at a point with H′J′ ≠ 0 is α₁² = α₂² = 1/(32π), V = 0 (**checked
symbolically**). So #161's statement that Theorem E's metrics are exact solutions with one normal and one phantom
scalar is now computed, and nothing in it changes. No discrepancy with #161, #65 or #78 was found.

## 1. Signs and conventions

Signature −+++. R_bc = ∂_aΓ^a_bc − ∂_cΓ^a_ab + Γ^a_adΓ^d_bc − Γ^a_cdΓ^d_ab, G_ab = R_ab − ½R g_ab. As in
`two_scalars.md`,

    T_ab = ∂_aφ₁∂_bφ₁ − ∂_aφ₂∂_bφ₂ − g_ab(½(∇φ₁)² − ½(∇φ₂)² + V).

A **normal** scalar enters T with +, a **phantom** with −. This is the issue's T₁ − T₂ with one potential, V = V₁ −
V₂. The code's stress function gives, for φ = ωt in Minkowski, T_tt = T_xx = +ω²/2 for a normal scalar and −ω²/2 for a
phantom (`conv_energy`), so a normal scalar has positive energy density. With a = √(8π)dφ₁ and b = √(8π)dφ₂, G = 8πT
is R_ab = 8πV g_ab + a_aa_b − b_ab_b, the form of `two_scalars.md`; `E_rform` checks it on the fields below, together
with a = ½(c dH − (4/c)dJ) and b = ½(c dH + (4/c)dJ), the identification of `planar_closure.md` §3.

## 2. The residual

For the metric of Theorem E (printed in `planar_residual_output.txt`):

| quantity | value |
|---|---|
| R_pq | −2H′J′ (from R = a a − b b, `E_rform`) |
| G_yy = G_zz | −2e^{−2HJ}H′J′ = −K, the other components 0 |
| ∂_pφ₁, ∂_qφ₁ | α c H′, −α(4/c)J′ |
| ∂_pφ₂, ∂_qφ₂ | α c H′, +α(4/c)J′ |

With α₁, α₂ free and V the value of the potential at a point, the nonzero components of G_ab − 8π(T₁ − T₂)_ab are, as
printed by `E_solve` (numerators, h = H′, j = J′, E = e^{2HJ}, A_i = α_i²):

| component | equation |
|---|---|
| pp | −8π c² h² (A₁ − A₂) = 0 |
| qq | −128π j² (A₁ − A₂) = 0 |
| pq | −8π E V = 0 |
| yy = zz | 2(16π(A₁ + A₂) h j + 4π E V − h j) = 0 |

Reading them: pp and qq give α₁² = α₂²; pq gives V = 0 (for a field of (p, q) alone (∂φ)² = 2φ_pφ_q/g_pq, so each
scalar has T_pq = 0, and G_pq = 0, which leaves only V in pq); yy then gives 32π α² = 1. The script confirms this
without reading: the 3 × 3 minors of the coefficient matrix, with the row denominators cleared, have a gcd 4096 E P³ h
j (P stands for π), and together with t h j E c P − 1 they generate the unit ideal (Gröbner basis [1]), so the matrix
has rank 3 at every point with h j E c ≠ 0. `linsolve` returns the single point (1/(32π), 1/(32π), 0), which makes
every equation vanish identically, and it does not contain c (`E_solve`).

With these values all 16 components of the residual simplify to exactly 0 for H(p), J(q) undefined and c a real
nonzero symbol (`E_residual`; G and T₁ are nonzero, so the zero is not vacuous). The same function with explicit
functions gives 0 for (H, J, c) = (p³ − p, q² + q, 3/2) and (sin p, q³ − q, −2/5), the second with H′J′ of both signs
(`E_explicit`).

## 3. For which c

For every c ≠ 0. The residual is a rational function of c, computed with c a free symbol, and it is 0, so it vanishes
wherever the fields are defined. c = 0 is not in the family (the fields have (4/c) J). For real c the fields are real.
Three remarks.

- c → −c sends φ₁ → −φ₁ and φ₂ → −φ₂, which leaves T unchanged, so c > 0 loses nothing.
- c is the field-space boost of `two_scalars.md` stage 1. With k = a + b = c dH and m = a − b = −(4/c) dJ, a boost by
  a constant χ sends k → e^χ k and m → e^{−χ} m, which is the same pair with c → e^χ c. a a − b b = −½(k m + m k) is
  unchanged (**proved**), so T₁ − T₂ has no c in it. The output shows it directly: 8π(T₁ − T₂) has two nonzero
  components, yy and zz, both −2e^{−2HJ}H′J′ = G, and d/dc of every component is 0.
- The identity holds for any sign of H′J′: α² is a positive constant, not tied to the sign of K = 2H′J′e^{−2HJ}. So
  the pair of scalars exists also where K changes sign, which matters for #170.

Assumed: H, J ∈ C², so that the metric is C². The residual contains only H′ and J′.

## 4. The potential

V = 0 along the solution. A constant V ≠ 0 does not work: it adds 8πV g_ab to the residual (control `potential`, V =
1/(80π): pq component −e^{2HJ}/10, yy = zz = 1/10). The map (p, q) → (φ₁, φ₂) has dφ₁ ∧ dφ₂ = 8α²H′J′ dp ∧ dq
(**proved**), so on the open set H′J′ ≠ 0 it covers an open set of field space, and V(φ₁, φ₂) = 0 on that set. Outside
the image V is unconstrained. The Klein–Gordon equations □φ_i = −2e^{−2HJ}∂_p∂_qφ_i = 0 hold since each φ_i is a sum
of a function of p and a function of q (**proved**; not run).

## 5. Checks and controls

Run of record: `run_final.run2` in the raw folder `I171-rpi5` of negative-mass-raw (stdout copied to
`planar_residual_output.txt`), exit 0: 8 of 8 checks pass in 11 s, then the 11 controls, 165 s in all, one process at
nice 10 on a Raspberry Pi 5, peak resident 64 MB (`/usr/bin/time -v` on the checks of the first version of the script,
raw `dev_checks.run2.time`). Run 1 (same script before the printing of the equations; same outcomes) is kept. From the
repository root, in the environment of [`setup.md`](../../setup.md): `python3 notes/warp/planar_residual.py`;
`--no-controls` for the checks alone, `--mutate=NAME` for one control.

| check | what it tests | result |
|---|---|---|
| `conv_energy` | normal T_tt = T_xx = +ω²/2, phantom −ω²/2 | pass |
| `planted_p1` | #78's P1, two normal scalars, V = −1/(8π), sheared coordinates: residual 0 | pass |
| `planted_p2` | #78's P2, two phantoms, V = +1/(8π): residual 0 | pass |
| `planted_flat` | Minkowski, φ₁ = φ₂ = t²x + sin y e^z, one normal and one phantom: residual 0, T₁ ≠ 0 | pass |
| `E_solve` | residual linear in (α₁², α₂², V); rank 3 wherever h j E c ≠ 0; unique solution (1/(32π), 1/(32π), 0), no c | pass |
| `E_residual` | α₁ = α₂ = 1/(4√(2π)), V = 0, H, J undefined, c symbolic: 16 components 0 | pass |
| `E_rform` | R_ab − (a_aa_b − b_ab_b) = 0 with λ = 0, and a, b equal ½(c dH ∓ (4/c) dJ) | pass |
| `E_explicit` | the same function with two explicit (H, J, c) | pass |

The planted checks run the same `residual` function as `E_residual`. Each control changes the code or an input once,
and its targets were computed in a development run (`dev_controls.run1.log`) and declared in the raw README before the
run of record.

| control `--mutate=` | change | checks that fail (all declared) | residual it produces (`E_residual`) |
|---|---|---|---|
| `scale11-normal` | φ₁ × 11/10 | `E_residual`, `E_rform`, `E_explicit` | pp −21c²H′²/400, qq −21J′²/(25c²), yy = zz 21e^{−2HJ}H′J′/100 |
| `scale11-phantom` | φ₂ × 11/10 | the same | pp +21c²H′²/400, qq +21J′²/(25c²), yy = zz 21e^{−2HJ}H′J′/100 |
| `flip-phantom` | T₁ + T₂ | `E_solve`, `E_residual`, `E_rform`, `E_explicit` | pp −c²H′²/2, qq −8J′²/c², yy = zz −2e^{−2HJ}H′J′ |
| `swap-roles` | φ₁ phantom, φ₂ normal | the same | yy = zz −4e^{−2HJ}H′J′ (= 2G_yy), pp = qq = 0 |
| `wrong-partner` | φ₂ = α(cH + (2/c)J) | the same | qq −3J′²/c², yy = zz −e^{−2HJ}H′J′/2 |
| `potential` | V = 1/(80π) | `E_residual` | pq −e^{2HJ}/10, yy = zz 1/10 |
| `ricci-sign` | Ricci tensor sign | `planted_p1`, `planted_p2`, `E_solve`, `E_residual`, `E_rform`, `E_explicit` | yy = zz +4e^{−2HJ}H′J′ |
| `conv-flip` | sign of T | `conv_energy`, `planted_p1`, `planted_p2`, `E_solve`, `E_residual`, `E_explicit` | yy = zz −4e^{−2HJ}H′J′ |
| `p1-sign` | one kinetic sign flipped in P1 | `planted_p1` | nonzero (first components printed) |
| `p2-scale11` | one P2 field × 11/10 | `planted_p2` | nonzero (printed) |
| `flat-flip` | phantom sign flipped in `planted_flat` | `planted_flat` | −8π(t⁴ + 4t²x² + e^{2z}) in tt (printed) |

All 11 fire on their declared targets, none with extras, none silent. The silences that were predicted and occur:
under `ricci-sign`, `planted_flat` (G = 0 in flat space) and `conv_energy` do not fail; under `conv-flip`,
`planted_flat` (T₁ − T₂ = 0 either way) and `E_rform` (it does not call the stress function) do not fail; the scalings
leave `E_solve` alone because its amplitudes are unknowns. The README's predictions (`dev_controls_predictions.txt`,
written before the development run) matched in all 11 cases.

The two controls the issue names, read from the printed residuals: scaling φ₁ by s = 11/10 gives −8π(s² − 1)T₁ and
scaling φ₂ gives +8π(s² − 1)T₂, with s² − 1 = 21/100 (pp entry ±21c²H′²/400 = ±8π · (21/100) · α²c²H′²); flipping the
phantom's sign gives −16πT₂ (pp entry −c²H′²/2). Both are nonzero wherever H′J′ ≠ 0 (yy = zz).

## 6. What `planar_closure.md` should now say

*Not edited here; for the coordinator.* §3, the paragraph after Theorem E's proof, last sentence ("*By hand; the
residual stage with these fields is not computed here.*"), could read:

> With α = 1/(4√(2π)) the scalars φ₁ = α(cH − (4/c)J) (normal) and φ₂ = α(cH + (4/c)J) (phantom), V = 0, satisfy G_ab
> − 8π(T₁ − T₂)_ab = 0 exactly for H, J undefined and every c ≠ 0 (*checked symbolically*,
> [`planar_residual.md`](planar_residual.md), #171). The normalizations and V = 0 are forced wherever H′J′ ≠ 0.

§8, "By hand: … the identification of the scalar fields" can drop the last item. The coordinator's review line "Not
computed: the residual stage …" is answered by #171.

## 7. What the code proves, and what is assumed

Exact (sympy, with undefined H(p), J(q) and a symbolic c, or exact explicit functions): the Ricci and Einstein tensors
of Theorem E's metric; the residual components for free α₁², α₂², V; the rank statement (Gröbner basis of the minors
with the saturation t h j E c P − 1); the single solution; the zero residual at the solution; R = λg + aa − bb with a,
b as in `planar_closure.md`; the same residual for two explicit functions; P1, P2 and the flat case through the same
function.

By hand: the dependence on c as a boost, the pq argument, the Jacobian of (p, q) → (φ₁, φ₂), the Klein–Gordon
equations.

Independence: the Einstein routine and the stress function were written for this script and share only sympy with
`phantom_source.py`, `two_scalars.py` and `planar_closure.py`, from which #161's K and the equivalence came. What they
share with the planted checks is the same functions, which is the point of the planted checks. Zero is tested by
`sympy.simplify(e) == 0`; sympy returning 0 is an exact vanishing, and no component needed any other test.

Assumed: two real, minimally coupled scalars of opposite sign with one potential V(φ₁, φ₂), G_ab = 8πT_ab, H, J ∈ C²,
sympy 1.14.0.

Not covered: Theorem E's global statements (they do not depend on the scalars). Other splittings of the same R: where
H′J′ ≠ 0, `two_scalars.md` stage 2 allows only a constant χ (ω = 0 there, since k and m are closed and independent),
so c is the only freedom; not rechecked here.

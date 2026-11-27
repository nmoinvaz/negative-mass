# Junction algebra of the S-brane and the Schwarzschild–de Sitter gluing (#224)

November 2026. Task T4 of `baby_universes.md` §4 (#211). Code `baby_universes_junction.py`, output `baby_universes_junction_output.txt`
(28 checks, 0 failed, 3 s). Raw `I224-rpi5` in negative-mass-raw: `README.md` (plan, by-hand values, declared mutation targets, reading
log), `run3*` (the run of record), `reading/` (the [FMM90] attempts). Status words as in AGENTS.md: *proved* (by hand) and
*checked symbolically* (sympy, inside the script). Nothing here says that a universe forms; each statement is about the algebra of a
gluing.

## Result

- **[BHR21] (4)–(5) hold** as the junction conditions between the Schwarzschild black-hole interior and a Kantowski–Sachs cosmology,
  with two readings fixed. The printed g is −f(r₀) = r_s/r₀ − 1, not the determinant that the text after (1) calls it (its square
  root has the dimension of length², not 1). For the action (1) at κ = 1 the right-hand side is S/2, a factor 2 below the printed S.
  The signs are as printed. *proved*, *checked symbolically* (P1–P3).
- **Positive tension is required** for the S² radius to expand, at every r₀ inside the horizon, and for the volume to expand where
  r₀ < 3 r_s/4. It is not required for ȧ/a, which is positive in the Schwarzschild past already (X1–X3).
- **NEC and WEC.** The S-brane (p = −σ, energy density 0) violates both for σ > 0, saturates both at σ = 0 and satisfies both for
  σ < 0. The vanishing energy density is forced by Einstein's equations. Any spacelike shell that reverses the contraction of the S²
  radius has p_x < 0 and violates the NEC (E1, E2).
- **Schwarzschild to de Sitter** (the setting of [FMM90]): the shell has p_x < 0 for every r₀ > 1/H, r_s and H; for r₀ < 0.9714 r_s both
  pressures are negative, so NEC and WEC fail in every direction. No shell-free matching exists. The isotropic S-brane glues the two
  only at H r₀ = √((F+1)(9F+1))/(3F+1), between 1 and 2/√3, and its tension is positive (D3–D5).
- **[FMM90]'s text was not reached** (section 8). Its abstract does not say what the layer's stress is, and this note does not
  claim to know.
- **Mutation.** Flipping the sign of the jump of K fires 8 of the 28 checks (J1 P0a P3 X3 E2 D3 D6 D6b); all 8 mutations fire exactly
  their declared targets.

## 1. Conventions

Signature (−+++), G_μν = κ T_μν, κ = 8πG ([BHR21] uses (+−−−); every physical statement below is signature independent). Hypersurface
r = r₀ with future-directed unit normal n; K^i_j = ½ h^{ik} L_n h_kj, so K^x_x = ȧ/a on a Kantowski–Sachs side whatever the signature.
[X] = X₊ − X₋, with + the future side. A shell has T^μ_ν = S^μ_ν δ(ℓ) with ℓ the proper time across it, S^i_j = diag(p_x, p_θ, p_θ) and
S^n_n = 0 (shown in section 2). An S-brane of tension σ, action −σ ∫√h d³x, has p_x = p_θ = −σ.

Inside the horizon f = 1 − r_s/r < 0. Write F := −f(r₀) = r_s/r₀ − 1 > 0. In the black-hole interior r is a time that decreases to the
future, n^r = −√F.

```
 future   de Sitter wedge, r increasing:   b = r = cosh(Hτ)/H,  a = sinh(Hτ)/H           (or a Kantowski–Sachs cosmology)
            │
 ═══════════╪═══════════  spacelike shell at r = r₀:  b continuous, [K^θ] > 0, p_x < 0
            │
 past     Schwarzschild interior, r decreasing:  b = r,  a = √F(r)
```

## 2. The junction condition (by hand)

Gaussian normal coordinates, ds² = −dτ² + A(τ)² dx² + B(τ)² dΩ², n = ∂_τ, K^x_x = Ȧ/A, K^θ_θ = Ḃ/B. The mixed Einstein tensor is (output)

```
G^τ_τ = −(2 Ȧ Ḃ/(A B) + Ḃ²/B² + 1/B²),   G^x_x = −(2 B̈/B + Ḃ²/B² + 1/B²),   G^θ_θ = −(B̈/B + Ä/A + Ȧ Ḃ/(A B)).
```

Only second derivatives can carry a δ(τ). G^τ_τ has none, so the shell has **no energy density** (S^n_n = 0), whatever it is made of.
With A and B continuous and Ȧ, Ḃ jumping, ∫Ä dτ = A [K^x] and ∫B̈ dτ = B [K^θ]:

```
κ p_x = −2 [K^θ],      κ p_θ = −([K^x] + [K^θ]),      i.e.   [K^i_j] − δ^i_j [K] = κ S^i_j.        (J)
```

This is Israel's thin-shell condition [Isr66] for a spacelike hypersurface; J1 below derives it without using it.
Inverse: [K^θ] = −κ p_x/2 and [K^x] = −κ p_θ + κ p_x/2. For an S-brane (p_x = p_θ = −σ) every K^i_j jumps by s = κσ/2.

Sign anchor: in flat FRW, Ḣ = −(κ/2)(ρ + p) gives [H] = −(κ/2) p_s for a layer with ∫p dτ = p_s; a negative-pressure layer speeds the
expansion up. (J) agrees (P0a). The Hamiltonian constraint, κρ = ½(R₃ + K² − K^i_j K^j_i) with R₃ = 2/r₀², gives across the shell
[ρ] = −S^i_j K̄^j_i, where K̄ is the average of the two sides.

Checks: J1 extracts the singular part of G^i_j from the kinked metric (no formula assumed); J2 applies (J) to K computed from the metrics;
P0a runs a flat FRW layer with the hand values above; D6 tests [ρ] = −S·K̄ on the de Sitter gluing.

## 3. The two sides

**Schwarzschild interior** (n^r = −√F). K^t_t = ½ n^r ∂_r g_tt/g_tt with g_tt = −f = F(r), ∂_r F = −r_s/r²; K^θ_θ = n^r/r:

```
K^x_x(−) = r_s/(2 r₀² √F) = (F+1)/(2 r₀ √F) > 0,      K^θ_θ(−) = −√F/r₀ < 0,       K(−) = (1 − 3F)/(2 r₀ √F).
```

The t direction stretches and the sphere shrinks; the volume r²√F is largest at r = 3r_s/4, where K(−) = 0. Computed twice (J2a: Lie
derivative in Schwarzschild coordinates; S0c: the interior rewritten as a Kantowski–Sachs cosmology with A = √F(r), B = r, dr/dτ = −√F,
which satisfies the vacuum equations). R₃ + K² − K^i_j K^j_i = 0 on this side (S0d).

**(i) Kantowski–Sachs future.** ds² = −dτ² + a² dx² + b² dΩ², K^x_x = ȧ/a, K^θ_θ = ḃ/b at τ₀. The induced metric is continuous when b(τ₀) = r₀
and a(τ₀)² dx² = F dt²; a rescaling of x moves a(τ₀) freely and leaves ȧ/a unchanged. (J) then fixes (ȧ/a, ḃ/b) from the shell stress.

**(ii) de Sitter future.** In the static form f₂ = 1 − H²r² < 0 for r > 1/H, r increases to the future, and the surfaces r = const are
R × S². In Gaussian normal form (the hyperboloid −X₀² + X₁² + … = 1/H² with (X₀, X₁) = c (cosh x, sinh x) and (X₂, X₃, X₄) = r n̂)

```
ds² = −dτ² + (sinh(Hτ)/H)² dx² + (cosh(Hτ)/H)² dΩ²,      G^μ_ν = −3H² δ^μ_ν  (S0b, D2b: both forms),   τ > 0 the expanding wedge.
```

The area radius r decreases in the Schwarzschild past and increases in the de Sitter future: the shell sits at a bounce of r. Matching
b(τ₀) = r₀ gives c := H r₀ = cosh(Hτ₀) ≥ 1, so r₀ ≥ 1/H, and with t := sinh(Hτ₀) = √(c²−1):

```
K^x_x(+) = H c/t = c²/(r₀ t),      K^θ_θ(+) = H t/c = t/r₀,      K^x K^θ = H² > 0.
```

The τ₀ < 0 wedge contracts and is not a new universe; it is not considered.

## 4. [BHR21] (4)–(5) against the derived forms

Printed: (4) ȧ/a + √g f′/(2f) = S, (5) ḃ/b + √g/r₀ = S, "where a prime denotes the derivative with respect to r".

| | printed | derived | |
|---|---|---|---|
| (4) | ȧ/a + √g f′/(2f) = S | ȧ/a − K^x(−) = s, with K^x(−) = r_s/(2r₀²√F); for g := −f(r₀) = F, √g f′/(2f) = −r_s/(2r₀²√F) | agree (P1) |
| (5) | ḃ/b + √g/r₀ = S | ḃ/b − K^θ(−) = s, K^θ(−) = −√F/r₀ | agree (P1) |
| g | "the determinant of the space-time metric" | −f(r₀) = r_s/r₀ − 1, dimensionless | notation (P2) |
| right side | S | s = κσ/2, which is S/2 for the action (1) at κ = 1 | factor 2 (P3) |

- g. With g read as the determinant (r⁴ sin²θ up to sign) √g has dimension length², and neither printed term scales as 1/length. The
  identification g := −f(r₀) makes both terms scale as 1/length and reproduces the derived K(−). The mutation `printed_g` fires P1 and P2.
- Factor 2. The action (1) has ½R, so κ = 1, and the brane term −S ∫√h gives p = −S; (J) then gives [K^i_j] = κσ/2 = S/2. The printed
  right-hand side is S. The sign and the structure of (4)–(5) are unaffected. The threshold for ḃ/b > 0 is S > √F/r₀ in the printed normalization and
  σ > 2√F/r₀ for the action (1) at κ = 1. The paper's own
  normalization of S in (4)–(5) is not stated there, and its details are in arXiv:2103.02842 (not read), so this is reported, not
  corrected.
- Orientation. The printed signs of the K(−) terms are those of the black-hole interior (r decreasing to the future). With the
  white-hole orientation (`orientation`) P1 fails.
- "Vanishing energy density" holds by (J) for any spacelike shell, not by a choice of the brane.

## 5. Which sign of the tension gives expansion

From (J) with n^r = −√F: ȧ/a = K^x(−) + s, ḃ/b = K^θ(−) + s, and the volume rate K(+) = K(−) + 3s. Thresholds (X1, X2; in σ = 2s/κ, X3):

```
ȧ/a > 0   iff  s > s_a = −(F+1)/(2 r₀ √F)   < 0              σ > −(F+1)/(κ r₀ √F)
ḃ/b > 0   iff  s > s_b =  √F/r₀              > 0              σ >  2√F/(κ r₀)
K(+) > 0  iff  s > s_v = (3F−1)/(6 r₀ √F)    > 0 iff F > 1/3, i.e. r₀ < 3 r_s/4
```

| r₀/r_s | F | r₀ s_a | r₀ s_b | r₀ s_v |
|---|---|---|---|---|
| 0.90 | 0.1111 | −1.6667 | 0.3333 | −0.3333 |
| 0.75 | 0.3333 | −1.1547 | 0.5774 | 0.0000 |
| 0.50 | 1 | −1.0000 | 1.0000 | 0.3333 |
| 0.10 | 9 | −1.6667 | 3.0000 | 1.4444 |
| 0.01 | 99 | −5.0252 | 9.9499 | 4.9582 |

So the statement "a positive S-brane tension is required in order to get expansion" is *proved* for the sphere radius at every r₀ < r_s,
and for the volume where r₀ < 3r_s/4, which includes the regime the paper has in mind, r₀ close to the singularity. It is not true of
ȧ/a: a negative tension down to s_a still leaves a expanding. "For a sufficiently large tension the future is expanding" holds: the
thresholds are finite. The post-brane energy density is κρ₊ = s (2K(−) + 3s) (D6b); since 2K(−) + 3 s_b = 1/(r₀ √F) > 0, it is positive
whenever s > s_b. A vacuum future (ρ₊ = 0) has s = 0 or s = −2K(−)/3 < s_b, so it cannot have ḃ/b > 0; the matter that the paper says the
brane's decay produces is required by the constraint, with this energy density.

## 6. Surface stress and the energy conditions

Orthonormal frame on the shell, null k = n + e with e a unit spatial vector, observer u = γ(n + v). For S^n_n = 0 and S^i_j = diag(p_x, p_θ, p_θ):

```
T(k,k) = p_x e_x² + p_θ (1 − e_x²),        T(u,u) = γ² (p_x v_x² + p_θ v_θ²),      so  NEC ⇔ WEC ⇔ p_x ≥ 0 and p_θ ≥ 0.
```

(WEC reduces to NEC here because the energy density is zero; E2.)

| σ | S^i_j | NEC, T(k,k) | WEC, T(u,u) | S² radius re-expands |
|---|---|---|---|---|
| > 0 | −σ δ^i_j | −σ < 0 for every null k: violated | −σ γ² v² ≤ 0, zero only for the comoving observer: violated | iff σ > 2√F/(κ r₀) |
| 0 | 0 | 0: saturated | 0: saturated | no |
| < 0 | −σ δ^i_j, positive | −σ > 0: satisfied | satisfied | no |

(E1.) The NEC statement is a statement about the coefficient of δ(τ). It confirms the paper's "violates the usual energy conditions",
for σ > 0, and shows that the violation is the price of the expansion, not a feature of the brane that could be traded away:

- By (J), [K^θ] = −κ p_x/2. A shell with p_x ≥ 0 (the NEC) lowers ḃ/b or leaves it. Since K^θ(−) = −√F/r₀ < 0, then ḃ/b(+) ≤ −√F/r₀ < 0:
  **no NEC-satisfying shell makes the S² radius re-expand**, isotropic or not. *proved*; E2.
- The same holds for any past and future with r decreasing then increasing: K^θ = n^r/r gives [K^θ] = (√|f₁| + √|f₂|)/r₀ > 0, so
  p_x = −2(√|f₁| + √|f₂|)/(κ r₀) < 0.

## 7. Schwarzschild to de Sitter

With (J), K(−) of section 3 and K(+) of (ii) (D2, D3), t = sinh(Hτ₀), c = H r₀ = √(1+t²):

```
κ r₀ p_x = −2 (t + √F)                              < 0   for all F, t > 0,
κ r₀ p_θ = −[ (2t²+1)/t + (F−1)/(2√F) ]             p_θ ≥ 0 only if F < (3−2√2)² = 17 − 12√2 ≈ 0.0294  (r₀ > 0.9714 r_s).
```

The second line: (2t²+1)/t − 2√2 = (√2 t − 1)²/t ≥ 0 and y² + 4√2 y − 1 = (y − (3−2√2))(y + 3 + 2√2) with y = √F (D3b). A sample, in units κ r₀ p
(output):

| F | r₀/r_s | t | H r₀ | κ r₀ p_x | κ r₀ p_θ | NEC, WEC |
|---|---|---|---|---|---|---|
| 0.0100 | 0.9901 | 0.7071 | 1.2247 | −1.6142 | +2.1216 | violated along x |
| 0.0294 | 0.9714 | 0.7071 | 1.2247 | −1.7573 | −0.0000 | violated along x |
| 1.0000 | 0.5000 | 0.7071 | 1.2247 | −3.4142 | −2.8284 | violated in all directions |
| 100.0000 | 0.0099 | 0.1000 | 1.0050 | −20.2000 | −15.1500 | violated in all directions |
| 100.0000 | 0.0099 | 0.7071 | 1.2247 | −21.4142 | −7.7784 | violated in all directions |

- **No shell-free matching** (D4, *proved*). On the Schwarzschild side K^x K^θ = −(F+1)/(2r₀²) < 0, on the de Sitter side K^x K^θ = H² > 0.
  [K^i_j] = 0 would need equal products. In the Schwarzschild interior the t direction stretches while the sphere shrinks; in the
  de Sitter wedge both expand together. A layer is needed for every r_s, r₀, H.
- **The S-brane of the paper is a one-parameter special case.** An isotropic tension needs [K^x] = [K^θ]. With t = sinh(Hτ₀) one has
  K^x(+) − K^θ(+) = 1/(r₀ t) and K^x(−) − K^θ(−) = (3F+1)/(2 r₀ √F), so isotropy holds iff t = 2√F/(3F+1), i.e. H r₀ = √((F+1)(9F+1))/(3F+1).
  Then s = 3√F (F+1)/(r₀ (3F+1)) > 0 (D5). H r₀ lies between 1 and 2/√3, the maximum at F = 1/3 (r₀ = 3r_s/4). For other H the shell
  has p_x ≠ p_θ.

| F | r₀/r_s | H r₀ | r₀ s | r₀ s_b |
|---|---|---|---|---|
| 0.01 | 0.9901 | 1.01868 | 0.2942 | 0.1000 |
| 1 | 0.5000 | 1.11803 | 1.5000 | 1.0000 |
| 100 | 0.0099 | 1.00221 | 10.0664 | 10.0000 |
| 10⁴ | 0.0001 | 1.00002 | 100.0067 | 100.0000 |

  The jump is s = s_b + t/r₀, above the threshold s_b of section 5 at every F, so the future has ḃ/b > 0 (it is de Sitter), and the S² radius at the
  shell is almost the de Sitter throat radius, H r₀ ≈ 1, for r₀ ≪ r_s.
- **Constraint** (D6): κρ(+) = ½(R₃ + K² − K^i_j K^j_i) = 3H² on the de Sitter side, 0 on the Schwarzschild side, and [ρ] = −S^i_j K̄^j_i holds
  for the derived S. This ties the shell stress to the vacuum energy of the future; it fails under a flipped jump or a dropped trace term.

What this does not settle: which r₀ and H a limiting-curvature argument selects, whether the future de Sitter wedge is the whole future
(the wedge has a(0) = 0, the bifurcation sphere), and anything about formation or stability.

## 8. [FMM90]

The text was not reached. INSPIRE (control number 267342) lists report number IC/88/91 (ICTP, Trieste) and a KEK scan id (KEKSCAN
2000-33-351), and OSTI has the abstract only. The KEK server answered "Request Rejected" to plain requests, the ICTP archive addresses
tried did not answer or gave 404, and CERN's server posed a bot challenge; no browser agent string was set and no shadow library was
tried. Every attempt is in raw `reading/fmm90_attempts.txt`. The paper needs to be supplied if the comparison is wanted; the KEK scan
(IC/88/91) is the likeliest free copy.

What the abstract says (OSTI and INSPIRE agree): the interior is studied "on the assumption that some limiting curvature exists"; "the
Schwarzschild metric inside the black hole can be attached to the de Sitter one at some spacelike junction surface which may represent
a short transition layer"; "the method of massive thin shells by Israel is used to obtain the characteristics of this layer"; "instead of
the singularity the closed world can be formed inside the black hole"; this "may also be valid in a more general case provided the
gravitation theory is asymptotically free and the limiting curvature exists"; the fate of the evaporating black hole is "briefly discussed".

What it does not say: the layer's stress, an energy condition, the matching parameters, or whether the de Sitter side is the R × S²
wedge used here. Section 7 is derived from the geometry alone and holds for any junction on r = r₀ to the expanding de Sitter wedge.
Whether [FMM90] obtain the same shell is not known.

## 9. Proposed text for the coordinator's files

`baby_universes.md` (the file is the coordinator's; nothing is edited here):

- §1 table, last cell of the [FMM89] [FMM90] [EB01] [BHR21] row: "junction algebra, done: #224".
- §2 row 7: "the junction algebra: positive tension is needed for the S² radius to expand and the shell violates the NEC and WEC | [BHR21],
  [FMM90] | done: #224, `baby_universes_junction.md` (printed (4)–(5) hold with √g = √(−f(r₀)) and a right-hand side S/2 for the
  action (1); any shell that reverses the contraction of the S² radius has p_x < 0)".
- §4, T4, after "Covered by": "**Done in #224.** (4)–(5) reproduced up to the notation of g and a factor 2 in the normalization of S;
  the NEC and WEC are violated for σ > 0; the Schwarzschild to de Sitter shell has p_x < 0. [FMM90]'s text was not available; the
  part on it is open until the coordinator supplies the paper."
- §6, "Not read": add "[BHR21]'s companion arXiv:2103.02842 settles the factor 2 and the symbol g, if it states them".

`literature.md`:

- [BHR21], after "Does not contain": "Junction check (#224, `baby_universes_junction.md`): (4)–(5) are the Israel conditions of the
  Schwarzschild interior and a Kantowski–Sachs cosmology with √g read as √(−f(r₀)); the printed right-hand side S is twice the one
  that the action (1) gives at κ = 1. Positive tension is needed for ḃ/b > 0 at every r₀ < r_s and for the volume at r₀ < 3r_s/4;
  the brane violates the NEC and WEC for σ > 0."
- [FMM89], [FMM90], after "Full text not read": "Reach (#224): INSPIRE gives IC/88/91 and KEK scan 2000-33-351; the KEK server rejected
  plain requests, the ICTP and CERN addresses tried failed, OSTI has the abstract only. A Schwarzschild interior glued to the
  expanding de Sitter wedge at r = r₀ needs a shell with p_x < 0 for every r₀, r_s, H (`baby_universes_junction.md` §7); whether
  [FMM90] state this is not known."
- References: [Isr66] W. Israel, "Singular hypersurfaces and thin shells in general relativity", Nuovo Cim. B 44S10, 1 (1966), erratum
  Nuovo Cim. B 48, 463 (1967), doi:10.1007/BF02710419, doi:10.1007/BF02712210 (INSPIRE record 50008 checked; the paper itself not read).

## 10. Checks, mutations and limits

Run: `python3 baby_universes_junction.py` from `notes/` (exit 1 if a check fails); `--mutate=all` runs the eight mutations below in
subprocesses and compares with the declared targets.

| id | what | status |
|---|---|---|
| S0a S0b S0a′ | Schwarzschild, de Sitter, flat FRW solve the equations (the Einstein routine, with its sign anchor) | pass |
| J1 J1b | singular part of G^i_j gives (J); G^τ_τ has none | pass |
| J2a J2b J2c S0c S0d | K(−) twice, K two ways, unit normal, constraint on the Schwarzschild side | pass |
| P0a P0b | planted: flat FRW layer (hand values); Schwarzschild to Schwarzschild, S·K̄ = 0 | pass |
| P1 P2 P3 | printed (4)–(5): equal with g := −f(r₀); dimensions; factor 2 | pass (P3 is the finding) |
| X1 X2 X3 | thresholds in s and in σ, and their signs | pass |
| E1 E2 | NEC and WEC of the S-brane and of a general shell; no re-expansion under the NEC | pass |
| D2 D2b D3 D3b D4 D5 D6 D6b | de Sitter K (two forms), shell stress, p_θ threshold, no shell-free matching, isotropic S-brane, constraint, post-brane density | pass |

| mutation | change | fires |
|---|---|---|
| `jump_sign` | [K] = K₋ − K₊ | J1 P0a P3 X3 E2 D3 D6 D6b |
| `trace_term` | the −δ^i_j [K] term dropped | J1 P0a P0b P3 X3 E2 D3 D6 D6b |
| `tension_sign` | p = +σ | P0a P3 X3 E1 |
| `orientation` | white-hole normal n^r = +√F | J2a S0c P1 X1 X2 X3 E2 D3 D5 |
| `lapse` | n^r = −F | J2c J2a S0c S0d P0b P1 X1 X2 X3 D3 D4 D5 D6b |
| `printed_g` | g read as −det g_μν | P1 P2 |
| `ds_swap` | de Sitter a and b exchanged | S0b D2 |
| `ds_orientation` | de Sitter normal reversed (static form) | D2b |

Independence: J1 (singular part of G) and J2 (the Israel formula) share only the metric ansatz and the Einstein routine, which S0a, S0b,
S0a′ test on known solutions. S0c and J2a are two computations of K(−). Every identity is exact (sympy); there is no tolerance.

Limits. A thick regularized layer was not integrated. The future of the shell was not evolved: nothing here says that a Kantowski–Sachs
future with these initial data isotropizes, which [BHR21] assume through an effective cosmological constant. The shell is a
distribution; "the NEC" is the statement on its δ coefficient. [FMM90] and arXiv:2103.02842 were not read.

References checked: [BHR21] arXiv:2105.07166 v1 (read, eqs (1)–(5) on the rendered page); [Isr66] INSPIRE record 50008 (bibliographic
data only); [FMM90] INSPIRE and OSTI abstracts only.

## Coordinator's review

Review of branch `wip/i224-junction` at 42a76a3 and raw `I224-rpi5` at 2517fef, recorded in negative-mass-raw `I224-mac-review`.

- **Scope.** The branch adds this note, the script and its output, nothing else.
- **Rerun.** On the Mac: 28 checks, 0 failed. The output is the same apart from the Python version and the timing.
- **Independent check.** The coordinator derived the Schwarzschild–de Sitter gluing by hand and coded it in sympy without
  this script (`ind224.py`).
  - Past side: the Schwarzschild interior, K^x = (F+1)/(2r₀√F) and K^θ = −√F/r₀.
  - Future side: de Sitter beyond its horizon, K^x = H²r₀/√G and K^θ = √G/r₀, with G = H²r₀² − 1.
  - Results:
    - κp_x = −2(√G + √F)/r₀ < 0 always, so no shell-free match exists;
    - the isotropic shell sits at H r₀ = √((F+1)(9F+1))/(3F+1), with maximum 2/√3 at F = 1/3;
    - its tension is positive.
  - These are section 7's results. The mutation that flips the de Sitter normal fires both checks.
- **Applied.** `baby_universes.md` (§1 table, §2 row 7, T4, §6) and `literature.md` ([BHR21], [FMM89]/[FMM90]) carry
  section 9's text, marked as from #224. The [Isr66] reference stays in this note only.
- **Open.** [FMM90]'s own text. The coordinator cannot reach it either, so the part on it stays open until the paper is
  supplied. [BHR21]'s companion arXiv:2103.02842 may settle the factor 2 and the symbol g. Neither is filed as a task.

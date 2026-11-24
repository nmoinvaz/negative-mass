# Independent review of #65's two-scalar result

November 2026, issue [#79](https://github.com/nmoinvaz/negative-mass/issues/79). The reviewed note is
[`two_scalars.md`](two_scalars.md). The reviewer had not seen `two_scalars.py`, its output, `phantom_source.py`,
`three_scalars.*` or `same_sign_scalars.*` when the checks below were written and run, and opened them only after the
final runs and their logs were committed to the raw repository. Section 9 says what the comparison showed.

What I recalled about #65 before starting: its verdict, and the claims of `two_scalars.md` read in this session. The
auto-memory index that loads automatically has one-line descriptions of #77 and #78 that mention a closure criterion and
a rank-1 condition for one-sign pairs. I read no file behind those lines, and no check below uses them.

Scripts written for the review, each with its stdout (all run from the repository root, in the environment of
[`setup.md`](../../setup.md); `--no-controls` runs the checks alone, `--mutate=NAME` one control; each exits 1 if a check
fails or a control does not fire):

| script | output | what it reviews |
|---|---|---|
| [`review_two_scalars_cartan.py`](review_two_scalars_cartan.py) | [`..._cartan_output.txt`](review_two_scalars_cartan_output.txt) | D for b(t, x, ρ) |
| [`review_two_scalars_lemma.py`](review_two_scalars_lemma.py) | [`..._lemma_output.txt`](review_two_scalars_lemma_output.txt) | the lemma |
| [`review_two_scalars_ideals.py`](review_two_scalars_ideals.py) | [`..._ideals_output.txt`](review_two_scalars_ideals_output.txt), [`..._ideals_certs_output.txt`](review_two_scalars_ideals_certs_output.txt) | coefficients, ideals, "no open set", accelerating bubble, b(ρ) |
| [`review_two_scalars_omega.py`](review_two_scalars_omega.py) | [`..._omega_output.txt`](review_two_scalars_omega_output.txt) | the ω construction, the control, b(ρ) closure |
| [`review_two_scalars_compare.py`](review_two_scalars_compare.py) | [`..._compare_output.txt`](review_two_scalars_compare_output.txt) | what `two_scalars_output.txt` prints, against the above (§9) |

Status words as in AGENTS.md. The raw folder is `I79-rpi5` in negative-mass-raw: its README holds the criteria, the
planted
cases and the control targets, all written before the first run, and every log, including the failed and the stopped ones.

## Verdicts

No error was found in #65's argument. Every check of #65's statements agrees. The gaps below are about wording and scope,
and none changes the verdict. Any change of status word is for the coordinator.

| item | #65's statement | verdict | status word |
|---|---|---|---|
| D | D = h(h² − P² − Q² + Y² + Z²) − P(Y² − Z²) + 2QYZ is det R₃ for b(t, x, ρ) | agrees (§1) | *checked symbolically*, two routes |
| lemma | R = λη + aa − bb iff λ = 0 and det R₃ = 0 | agrees (§2) | *proved* (re-derived by hand), 80 random tensors and 400 (R, λ ≠ 0) pairs |
| coefficients | 16r³D is a cubic in c², with the printed c⁰ and c⁶ | agrees (§3) | *checked symbolically* |
| ideal | the four coefficients generate ideal(2w − rw′, (w′ + rw″)²) when w′r ≠ 0 | agrees, with the localization said (§3) | *checked symbolically*, explicit cofactors |
| no open set | the two relations hold on no interval with w′ ≠ 0 | agrees (§3) | *proved*, two proofs |
| accelerating bubble | R₃ depends on w̃, its derivatives and A; with A w̃′ r ≠ 0 the seven coefficients generate (1) | agrees (§4) | *checked symbolically*, explicit cofactors and a hand chain |
| b(ρ) | D = h(h² − Q²); two solution families; both fail closure | agrees (§5) | *proved* by hand, *checked symbolically* |
| ω | ω = dm(K, M) k + ι_K dk is the only candidate | agrees (§6) | *proved* by hand, *checked symbolically* on planted data |
| control | Ellis–Bronnikov with a normal scalar passes all four stages | agrees on the algebraic, closure and residual stages (§6) | *checked symbolically* |

## 1. D from the Cartan formulas

#59's [`review.md`](review.md) §1 gives, for any b(t, x, y, z) with unit lapse and flat slices and n = ∂_t − b∂_x,

    R_00 = n(b_x) − b_x² − ½|∇⊥b|²,   R_xx = −n(b_x) + b_x² − ½|∇⊥b|²,   R_AB = ½ b_A b_B,
    R_0x = −½Δ⊥b,   R_0A = ½ ∂_x b_A,   R_xA = −½ n(b_A) + b_x b_A.

For b = b(t, x, ρ) in the frame (n, ∂_x, e_ρ), with ρ-derivatives along e_ρ and b_φ = 0, this is R₃ = [[P − h, Q, Z],
[Q, −P − h, Y], [Z, Y, h]], with P, Q, Z, Y, h as in #65. Its determinant, expanded by hand,
is (P − h)(−Ph − h² − Y²) − Q(Qh − YZ) + Z(QY + (P + h)Z) = h(h² − P² − Q² + Y² + Z²) − P(Y² − Z²) + 2QYZ. This is #65's
D.

`review_two_scalars_cartan.py` does not rely on those formulas:

- `blocks_jets`, `D_cartan_jets`: the imported Cartan routine of `review_cartan.py` (Jets, `axisymmetric_shift`,
  `cartan_ricci`), for b = G(t, x, y² + z²) with every jet a free symbol. P, Q, Z, Y, h are built here by total
  derivatives.
  Every entry of R₃ equals the block, R(e_φ, ·) = 0, and det R₃ equals D(printed), exactly (39 terms in the jets of G).
- `D_coordinate`: my own coordinate-Christoffel Ricci at 9 exact rational points of three explicit shifts b(t, x, ρ).
  Every
  entry of R₃ and its determinant equal the block and D, and the frame is orthonormal.
- `planted_D`: b = kρ² gives D = h(h² − 4k²) by both routes, and the four b(ρ) solutions of §5 give D = 0.
- Both Ricci routines first reproduce de Sitter and Ellis–Bronnikov (`cartan_planted`, `coord_planted`).
- `axis_D`: on the axis R(E_2, ·) = R(E_3, ·) = 0, R₀₀ + R₁₁ = 0 and det R₃ = 0. #65's "R₃ has a zero row, so det R₃ = 0
  and the
  form holds with λ = 0" also needs R₀₀ + R₁₁ = 0 there, so that R₃ is indefinite (the lemma's weighted trace). It
  holds, and
  `phantom_source.md` checks it as `identity_axis`; the note could cite it.

**A correction made during the review (claim A1 is false).** My first account of why a control (`nonaxi`) was silent was
that
"det R(n, ∂_x, ∂_y) at z = 0 equals the printed D for every b(t, x, y, z)". That claim is **false**: my check
`D_slice_nonaxi` failed 0 of 9 (`cartan.run2.log`, kept). For a general b, R₀₀ = P − h − b_z²/2 and R₁₁ = −P − h −
b_z²/2, with
h = b_y²/2 (`D_slice_offslice`: at 9 points with z ≠ 0 only these two entries of the six differ from the block, each by
exactly
−b_z²/2, and det R(n, ∂_x, ∂_y) ≠ D(printed)). The block has the printed form where b_z = 0 (`D_slice_bz0`: 9 of 9 for
shifts
even in z, which are not functions of ρ, at z = 0). The control `nonaxi` (b = G(t, x, y² + 2z²)) was silent on
`D_cartan_jets`
because b_z = 4zG_001 vanishes on the slice z = 0 that check uses. **For #65:** an axisymmetric b(t, x, ρ) has b_z = 0
on z = 0
exactly in the jets (`D_slice_bz0`, last clause), so #65's D is unaffected (*checked symbolically*). The coordinator
approved
replacing `D_slice_nonaxi` by the two checks after it failed.

## 2. The lemma

The statement: R real symmetric on an orthonormal frame, R_a3 = 0 for every a (including R₃₃ = 0) and R₀₀ + R₁₁ + 2R₂₂ =
0.
Then R = λη + aa − bb, λ real, a and b real covectors, iff λ = 0 and det R₃ = 0.

**Proof by hand (re-derived; it agrees with #65's).** Write M = R − λη = N ⊕ (−λ), N = R₃ − λ diag(−1, 1, 1), so that
N₀₀ + N₁₁ + 2N₂₂ = −2λ under the hypotheses.

*Sublemma.* The form q = a² − b² is not definite on any subspace U of dimension 2. If q < 0 on U∖0 then |b(u)| > |a(u)|
≥ 0, so b|_U is
injective into ℝ, which is impossible for dim U = 2. If q > 0, a|_U is injective, the same.

- λ > 0. M(e₃, e₃) = −λ < 0. If N(v, v) < 0 for some v in span(e₀, e₁, e₂), then M(se₃ + tv, ·) = −λs² + t²N(v, v) is
  negative
  definite on span(e₃, v), contradicting the sublemma. So N ⪰ 0, its diagonal is ≥ 0, and N₀₀ + N₁₁ + 2N₂₂ ≥ 0. That
  contradicts −2λ < 0.
- λ < 0. The same with the signs reversed.
- λ = 0. M = R₃ ⊕ 0 and the weighted trace R₀₀ + R₁₁ + 2R₂₂ vanishes. If R₃ ⪰ 0 or ⪯ 0 its diagonal vanishes, so R₃ = 0.
  Otherwise R₃
  has eigenvalues of both signs. If the form holds, rank M ≤ 2, so det R₃ = 0. Conversely, det R₃ = 0 with R₃ ≠ 0 forces
  the
  eigenvalues (μ₁ > 0, μ₂ < 0, 0), and a = √μ₁ e₁, b = √|μ₂| e₂ give R₃ = aaᵀ − bbᵀ. For R₃ = 0 take a = b = 0. ∎

The proof uses only that the kinetic form has signature (1, 1), so it covers any two-field sigma model with a Lorentzian
target metric,
as #65 says. For two fields of one sign, M ⪰ 0 and the first branch changes: λ < 0 with N ⪰ 0 of rank ≤ 1 is not excluded.

`review_two_scalars_lemma.py` tests the statement without the proof's route. It takes the real roots of the quartic
det(R − λη)
(not the gcd of the 3 × 3 minors) and counts eigenvalue signs of R − λη at 80 digits at each root.

- `ricci_from_T`: R_ab = 8π(∂φ₁∂φ₁ − ∂φ₂∂φ₂) + 8πVη_ab from G = 8πT, in symbols.
- `weighted_trace`: N₀₀ + N₁₁ + 2N₂₂ = −2λ.
- `search_planted`: 30 planted λ₀η + aa − bb (λ₀ = 0 in a third), in random rational Lorentz frames, are recovered;
  planted (+,+), (−,−)
  and (+,−,+) forms and 20 generic symmetric 4 × 4 give the empty set.
- `lambda_branch`: 400 random (R, λ ≠ 0) pairs satisfying the hypotheses never have inertia within (1, 1).
- `lemma_if`: 40 tensors with det R₃ = 0 (38 with R₃ ≠ 0) have λ-set {0}, and a, b rebuild R to 3 × 10⁻⁷⁹.
- `lemma_onlyif`: 40 tensors with det R₃ ≠ 0 have the empty λ-set.

One wording point. #65 says each tensor is tested "in a random rational Lorentz frame". The hypotheses are
frame-specific, so the
random frame tests only that the λ-search is Lorentz-covariant, which it is. It does not test the lemma in a second frame.

## 3. Constant speed: coefficients, ideals, no open set

b = w(r) − v, w = v(1 − f), c = cos θ. `r3_const`: R₃ from my own jet code (differentiating b(ξ, y, z) with r an
independent
symbol) equals the note's table and the imported Cartan Ricci of b = −vH((x − vt)² + y² + z²), exactly.

**Coefficients (`coeffs_const`).** 16r³D is a cubic in c² with no odd powers. With u = c², 16r³D = c₀ + c₁u + c₂u² +
c₃u³ and

    c₀ = 2 r w′² (r² w′⁴ − 4 w² w′² − (w′ + r w″)²)                          equal to the printed c⁰,
    c₃ = −2 (r w w″ + r w′² − w w′)² (2 r w w″ + r w′² − 2 w w′)              equal to the printed c⁶.

c₁ and c₂ are printed in the output (§2 there). With A = 2w − rw′ and B = w′ + rw″ the first is
c₀ = −2 r w′² (w′² A² + 2 r w′³ A + B²).

**The ideal (`ideal_const`).** In ℚ[r, w, w′, w″][1/(r w′)], ideal(c₀, c₁, c₂, c₃) = ideal(A, B²).

- Each c_k lies in ideal(A, B²): after the substitution w = (A + r w′)/2, w″ = (B − w′)/r no c_k has a monomial 1 or B,
  so c_k = α_kA + β_kB²
  with the cofactors printed by the script.
- A and B² lie in ideal(c_k). Explicit cofactors of degree 5 in (A, B), 39 nonzero coefficients each, were found by
  linear algebra over
  ℚ(r, w′) (no Gröbner basis) and verified by exact expansion. Every denominator is a monomial in r, w′, so the equality
  holds after
  inverting r w′ and not only over the field. The cofactors are in `review_two_scalars_ideals_certs_output.txt`.

So the square in B² is right, and the statement needs the localization, which #65 has as "with w′r ≠ 0". Over ℚ[r, w,
w′, w″] itself it is
not claimed and not tested. For a real point with r w′ ≠ 0, all four coefficients vanish iff A = 0 and B = 0
(`zero_set_const`: planted jets
with A = B = 0 give D = 0 at 10 random c, and A = 0, B ≠ 0 or B = 0, A ≠ 0 give D ≠ 0).

**No open set (`no_open_set`).** #65: differentiating 2w = rw′ gives w′ = rw″, which with rw″ = −w′ gives w′ = 0.
Checked: dA/dr = 2w′ − B
exactly. A second proof, by solving: A = 0 on an interval gives w = Cr², and B = 0, that is (rw′)′ = 0, gives w = K ln r
+ L. The functions
r², 1, ln r have Wronskian 4/r ≠ 0, so Cr² = K ln r + L on an interval forces C = K = L = 0, and then w′ = 0. Both
proofs agree.
*Proved.* The step from the zero set having no interior in (r, θ) to "no source on an open set meeting the wall" is as
in #65, and I find
no gap in it: a box in (r, θ) is open in the wall, D does not depend on t or ϕ, a cubic in c² vanishing on a θ-interval
is zero, and v ≠ 0.

## 4. The accelerating bubble

For b = −V(t)H((x − X(t))² + y² + z²), my comoving derivation uses n = ∂_t|_ξ − w̃∂_ξ with w̃ = X′ − V f and ∂_t|_ξ g =
A g for every spatial
derivative g of b, A = V′/V. Then P gains +A w̃′ c and Y gains −A w̃′ s/2, and Q, Z, h do not change. `acc_R3`: this
equals the imported Cartan
Ricci in lab coordinates, for generic jets of V, X, H, exactly in every entry. `acc_coeffs`: 16r³D has degree 6 in c,
the odd coefficients are
divisible by 4 A w̃′ r, and at A = 0 the seven reduce to §3's four.

**The unit ideal (`acc_unit`).** (A is the acceleration V′/V here, not §3's 2w − rw′.) With A w̃′ r inverted the seven
coefficients generate (1). Certificate: cofactors of degree 3 in (w, ω),
ω = (w′ − r w″)/(r w′²), 25 nonzero coefficients on the support c₀, c₄, c₅, c₆ with denominators monomials in r, w′, A,
verified by exact expansion
also in (w, w″). A hand chain that uses neither a Gröbner basis nor a search, each step exact in the output:

1. With Λ = r w w″ + r w′² − w w′ and ρ₁ = r w′²: c₆ = −2Λ²(2Λ − ρ₁) and c₅/(4Ar w′) = Λ(3Λ − ρ₁). Hence
   c₆ = (2ρ₁²/9)Λ + (c₅/(4Ar w′))(2ρ₁/9 − 4Λ/3), and Λ lies in the ideal because ρ₁ is a unit.
2. Modulo Λ, w(w′ − r w″) = ρ₁ − Λ, so w is a unit, and Λ = ρ₁(1 − wω).
3. Modulo Λ, c₄ becomes 2r²w′⁴F/w² and c₀ becomes 2r w′³G₀/w², with F a cubic and G₀ a quartic in w. Their resultant is
   Res_w(F, G₀) = 512 A² r⁴ w′⁹, a unit, so ideal(F, G₀) = (1).

The `acc-drop-odd` control drops c₅ and finds no certificate up to degree 5.

## 5. A shift b(ρ)

For b depending on ρ only, P = Z = Y = 0 by differentiation (`brho_D`), and Δ⊥b = b″ + b′/ρ, so R₃ = [[−h, Q, 0], [Q,
−h, 0], [0, 0, h]] and
D = h(h² − Q²). With u = b′ the equation b″ + b′/ρ = εb′² is the Bernoulli equation u′ + u/ρ = εu², whose solutions are
u = 1/(ρ(C − ε ln ρ)), so b = −ε ln|C − ε ln ρ| + b₀. Both signs of ε and C = 0, 3 solve it and have D = 0, exactly
(`brho_D`).
The classification up to x → −x, x → x + b₀t and a homothety is as in #65. On those solutions Q = −εh, so R = h((dρ)² −
((1 + εb)dt + ε dx)²).

## 6. The ω construction and the control

k = a + b, m = a − b. A field-space boost multiplies them by e^χ and e^{−χ}. Closure is ω = dχ with ω ∧ k = −dk, ω ∧ m =
dm, dω = 0.

**Proof by hand of the formula.** Take K, M with k(K) = m(M) = 1 and k(M) = m(K) = 0. Contracting ω ∧ k = −dk with K gives
ω(K)k − ω = −ι_K dk, so ω = ω(K)k + ι_K dk. Applying that to M and using k(M) = 0 gives ω(M) = dk(K, M). Contracting ω ∧
m = dm with M and
applying to K, with m(K) = 0, gives ω(K) = dm(K, M). So ω = dm(K, M) k + ι_K dk, unique when k ∧ m ≠ 0. ∎ A solution
exists iff the candidate
satisfies both 2-form equations and dω = 0. The choice of K, M is free modulo ker k ∩ ker m, and the candidate does not
depend on it when a solution exists.

- `omega_planted`: 8 planted pairs k = e^{−χ}(dφ₁ + dφ₂), m = e^{χ}(dφ₁ − dφ₂) with random polynomial φ₁, φ₂, χ = ln 2 +
  polynomial. The candidate
  equals dχ at the point for two different (K, M), solves both equations, and dω = 0 for the symbolic field on two of
  them.
- `omega_planted_fail`: 4 pairs with k ∧ dk ≠ 0 have no ω, and the candidate violates an equation.
- `omega_eb`: Ellis–Bronnikov, r² = (1 − β cos 2κl)/(2κ²), ultrastatic. My coordinate Ricci gives R_ab = λg_ab + aa − bb
  with λ = 2κ², a = √2κ dt and
  b = √(8π)φ₂′ dl, an identity in κ, β, l (cos 2κl and sin 2κl reduced with S² = 1 − C²). φ₂′ = C/r² from the arctan
  form to 80 digits at 6 points.
  G − 8πT = 0 symbolically. ω = 0 for the closed split. After a planted boost χ = tl + l²/3, ω = −dχ at 3 points to 50
  digits.
- `omega_brho`: for b″ + b′/ρ = εb′², exact in (b, b′, ρ): R_ab = h((dρ)² − ((1 + εb)dt + ε dx)²), the candidate
  violates an equation for both signs of √h,
  and (ν ∧ dν)_{txρ} = b′ for ν = dρ + (1 + εb)dt + ε dx. So no e^χ ν is closed (Frobenius).

#65's potential-stage argument is also correct as stated: for closed a and b, ∇^aG_ab = 0 gives ∇_bλ = (∇·a)a_b −
(∇·b)b_b, so dλ ∧ a ∧ b = 0.
I re-derived it by hand and did not put it in code.

## 7. Controls

27 controls in the four review scripts and 5 in the comparison script. In their final runs all 32 fire on their listed
targets. One was silent on a
named target in its first run:

- **`nonaxi` was silent on its named target `D_cartan_jets`; it fired `blocks_jets` and `planted_D`; the target list was
  changed after the run.**
  The reason is in §1, where the first explanation I gave is recorded as false. `D_cartan_jets` keeps controls that fire
  on it:
  `d-qyz-sign`, `d-drop-z2`, `blocks-q-sign`, `cartan-lapse` and `cartan-sign` (logs `cartan.mut_*.run2.log` and the
  final output).
- `acc-drop-odd` ran 12 min in a symbolic solve with no solution and was stopped by me. The control now first solves the
  same system at a
  numeric point; a certificate with monomial denominators would specialize, so a failure there is a failure everywhere.
  Its target did not
  change (`acc_unit`). The partial log, the change and the reason are in the raw README (A2). The same function serves
  the real check.
- Two of my own checks failed on a first run for bugs in the checks (`brho_D`: an abstract function with sympy's Abs; a
  format string in the
  lemma script; a keyword-argument bug in the ω script). They were fixed and rerun with the logs kept (raw README, A3).
  No criterion was changed.

| script | checks | controls | checks alone | whole run | peak memory |
|---|---|---|---|---|---|
| `review_two_scalars_cartan.py` | 10 | 8 of 8 fired | 12 s | 2 min 7 s | 71 MB |
| `review_two_scalars_lemma.py` | 6 | 5 of 5 | 17 s | 1 min 52 s | 72 MB |
| `review_two_scalars_ideals.py` | 9 | 8 of 8 | 19 s | 2 min 55 s | 72 MB |
| `review_two_scalars_omega.py` | 4 | 6 of 6 | 8 min 5 s | 10 min 31 s | 92 MB |
| `review_two_scalars_compare.py` | 6 | 5 of 5 | 2 s | 12 s | 61 MB |

Times are on a Raspberry Pi 5 (aarch64, 4 cores shared with other jobs), one process at a time at nice 10. Memory is the
peak resident size of the
whole run, from `/usr/bin/time -v`.

## 8. Gaps

None changes a verdict. They are reported, not fixed, in the note under review.

1. **The ideal needs its localization.** "The four coefficients generate the same ideal as 2w − rw′ and (w′ + rw″)²"
   holds in
   ℚ[r, w, w′, w″][1/(r w′)] (§3). The note says "with w′ r ≠ 0", which is the right condition. It could say where the
   ideal lives.
2. **The axis needs the second identity.** "On the axis R₃ has a zero row, so det R₃ = 0 and the form holds with λ = 0"
   also uses
   R₀₀ + R₁₁ = 0 on the axis (§1), so that R₃ is indefinite or zero.
3. **The random Lorentz frame in `lemma_samples`** tests the covariance of the λ-search, not the lemma in another frame
   (§2).
4. **The control's signs hold on all of 0 < β < 1.** "Numerical evidence only" in #65 covers the signs of the splitting,
   sampled on
   0 < β < 1/2. The split is an identity: R_ab = λg_ab + aa − bb with a = √2 κ dt and b_l² = (1 − β²)/(2κ²r⁴) > 0, exact
   in κ, β and l (§6). So the
   sign of each square follows by hand for every 0 < β < 1, and that sentence could go.
5. **`closure` divides by k·k m·m − (k·m)²** (no guard). It is zero only if k ∥ m with both nonzero, which `split_form`
   never returns (a rank-one
   form gets a zero covector, which `closure` treats separately and correctly: if cosh χ a and sinh χ a are closed then
   da = 0). So it is unreachable.
6. **A wrong first account of mine, now corrected (§1).** A claim about D at z = 0 for non-axisymmetric b was false and
   is on record. It does not
   concern #65, whose b is axisymmetric.

## 9. What the comparison with two_scalars.py showed

Made after the final runs and their commit, by reading `two_scalars_output.txt` and `two_scalars.py`, and by
`review_two_scalars_compare.py` (6 of 6 checks pass, 5 of 5 controls fire). Everything printed agrees with the review.

- The printed c⁰, c², c⁴, c⁶ of 16r³D equal mine exactly, and so do the printed odd coefficients c¹, c³, c⁵ of the
  accelerating bubble.
- At the sample points of its section 5, my D is nonzero at the 16 off-axis points and zero at the 7 on-axis ones. The
  tanh and bump points are
  recomputed from the profiles at 50 digits. The 4 random-jet points use its seed-57 draws in its order (read from its
  code, so this part is not
  independent of it). The printed 16 ± 4i and 4 ± 2i on the axis follow from my table (R₀₀ = −(ww″ + w′²), R₀₁ = −w′/r
  at c = 1).
- The printed b_l and λ of the control equal mine at 6 random parameter sets. The printed R₃ for b(ρ) equals mine.
- **Same formula, different code.** Its `closure` computes K and M from the Gram inverse of (k, m) and the candidate
  dm(K, M)k + ι_K dk, the formula
  of §6. It was written from the same derivation as mine. It also treats a single square separately.

Method differences, which make the two reviews complementary:

| step | two_scalars.py | this review |
|---|---|---|
| R₃ and D | #57's spherical and lab routes (`phantom_source.py`) | Cartan jets (`review_cartan.py`) and my coordinate Christoffel route |
| lemma | 24 tensors through its algebraic stage (gcd of the 3 × 3 minors) | 80 tensors and 400 (R, λ) pairs through the quartic det(R − λη) and eigenvalue signs |
| ideals | grevlex Gröbner bases with a saturation variable, membership both ways | explicit cofactors verified by expansion, and a resultant |
| "no open set" | elimination of W₀, W₂ and differentiation | the two ODE solution families and a Wronskian, and the same differentiation |
| accelerating bubble | route D of #61, Gröbner basis {1} | my comoving derivation against Cartan in lab coordinates, certificate and a hand chain |
| b(ρ) | the two solutions with C = 0, closure by the algebraic stage | every C, exact in (b, b′, ρ), closure by the candidate and by Frobenius |
| ω | `closure` | the formula by hand, planted data, the control, b(ρ) |

Not reviewed here: the algebraic-stage code itself (`candidates`, `split_form`, `sign_of`, whose sign decisions use
sampling where sympy cannot decide),
the exact-point tests beyond those listed above, the controls `eb2-ads`, `eb2-deform`, `eb2-unclosed`, `eb2-potential`,
`eb2-residual`, `planted-*`,
`acc-static` and `alc-*`, the potential stage in code, and the prior-work paragraph.

## 10. The lemmas for Lean (#135)

Statements only. `η = diag(−1, 1, 1, 1)`, `R₃` the block on indices 0, 1, 2. Each is stated over `ℝ`.

1. **Trace reversal** (`ring`). For `p q : Fin 4 → ℝ`, `V : ℝ` and `T_ab = p_a p_b − q_a q_b − η_ab (½ pᵀη⁻¹p − ½ qᵀη⁻¹q
   + V)`:
   `G = 8πT` implies `R_ab = G_ab − ½ (tr_η G) η_ab = 8π(p_a p_b − q_a q_b) + 8πV η_ab`.
2. **Determinant** (`Matrix.det_fin_three`, `ring`). `det !![P − h, Q, Z; Q, −P − h, Y; Z, Y, h] = h(h² − P² − Q² + Y² +
   Z²) − P(Y² − Z²) + 2QYZ`.
3. **Sublemma.** For linear `a b : V →ₗ[ℝ] ℝ` and a subspace `U` of dimension 2, `¬(∀ u ∈ U, u ≠ 0 → a u² − b u² < 0)`
   and the same with `> 0`.
   Proof: `b|_U` (resp. `a|_U`) is injective into `ℝ`.
4. **Lemma.** For `R : Matrix (Fin 4) (Fin 4) ℝ` with `R.IsSymm`, `∀ i, R i 3 = 0` and `R 0 0 + R 1 1 + 2 * R 2 2 = 0`:
   (a) if `∀ i j, R i j = lam * η i j + a i * a j − b i * b j` then `lam = 0` and `det R₃ = 0`;
   (b) if `det R₃ = 0` then `∃ a b, ∀ i j, R i j = a i * a j − b i * b j`.
   Part (a) uses the sublemma and the weighted trace `N₀₀ + N₁₁ + 2N₂₂ = −2λ`. Part (b) needs the spectral theorem for a
   symmetric 3 × 3 matrix
   whose weighted trace vanishes, which makes it indefinite or zero.
5. **Ricci blocks (assumed, not proved here).** For `b(t, x, ρ)` with unit lapse and flat slices, the frame Ricci block
   is the matrix of statement 2
   with `P, Q, Z, Y, h` as in §1, and `R(e_φ, ·) = 0`. This is the interface to differential geometry: it is a hand
   derivation (#59) with the
   Cartan and coordinate checks of §1 as evidence. Name it as a hypothesis.
6. **Coefficients** (`ring` after `s² = 1 − c²`). `16 r³ D` equals `Σ_k c_k (c²)^k` with `c₀, c₃` as in §3 and `c₁, c₂`
   as printed in the output.
7. **Ideal, constant speed.** For `r ≠ 0`, `w′ ≠ 0`: `c₀ = c₁ = c₂ = c₃ = 0 ↔ 2w = r w′ ∧ w′ + r w″ = 0`. (←) each `c_k
   = α_k A + β_k B²`;
   (→) `A = Σ γ_k c_k` and `B² = Σ δ_k c_k`, by `field_simp; linear_combination` with the 39-term cofactors of
   `review_two_scalars_ideals_certs_output.txt`, then `pow_eq_zero`.
8. **No open set.** For `I ⊆ Ioi 0` an open interval, `w : ℝ → ℝ` with `ContDiffOn ℝ 2 w I`, and `∀ r ∈ I, 2 w r = r *
   deriv w r ∧ deriv w r + r * deriv (deriv w) r = 0`:
   `∀ r ∈ I, deriv w r = 0`. Proof: differentiate the first relation and add the second.
9. **Polynomial step.** A real polynomial of degree ≤ 6 that vanishes at infinitely many values of `c²` (or `c`) is zero
   (`Polynomial.eq_zero_of_infinite_isRoot`).
10. **Accelerating bubble.** For `r ≠ 0`, `w′ ≠ 0`, `A ≠ 0`, the seven coefficients do not all vanish: `1 = Σ g_k c_k`
    with the 25-term cofactors of
    the certificates output, on the support `c₀, c₄, c₅, c₆`, after `field_simp`.
11. **b(ρ).** If `u′ + u/ρ = ε u²` on `(0, ∞)` with `u ≠ 0` then `u = 1/(ρ(C − ε ln ρ))` for a constant `C` locally, so
    `Q = −εh` and `D = 0`.
    For `ν = dρ + (1 + εb)dt + ε dx`, `(ν ∧ dν)_{txρ} = b′`. If `d(fν) = 0` then `0 = ν ∧ d(fν) = f ν ∧ dν`, so `f = 0` where `b′ ≠ 0`: no integrating factor.
12. **The ω formula** (pure linear algebra on `Fin 4`, `Finset.sum`). For `k m K M ω : Fin 4 → ℝ`, antisymmetric `dk
    dm`, `k·K = m·M = 1`, `k·M = m·K = 0`, if
    `∀ i j, ω i * k j − ω j * k i = −dk i j` and `∀ i j, ω i * m j − ω j * m i = dm i j`, then `ω j = (Σ K i M l dm i l) k j + Σ K i dk i j`.
13. **Potential stage** (by hand here). For closed `a, b` and `R = λg + aa − bb`, `∇^a G_ab = 0` gives `∇_b λ = (∇·a)a_b
    − (∇·b)b_b`.

## 11. What two_scalars.md should say about its review status

Suggested sentence: "Reviewed independently in #79 ([`review_two_scalars.md`](review_two_scalars.md)): D, the lemma, the
constant-speed coefficients and ideals, the
accelerating bubble, b(ρ) and the ω construction were re-derived or re-checked with code that shares only sympy and
`review_cartan.py`'s Cartan
routine with `phantom_source.py`. No error was found. The lemma is *proved* by hand and was re-derived; the rest is
*checked symbolically*."
The status words in `two_scalars.md` stand. The wording points of §8 (1, 2, 4) are optional. The sampled-sign sentence
for the control can be dropped.

## Coordinator's review

Raw `I79-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| scope | the twelve files the issue allows, and raw `I79-rpi5` |
| rerun the five scripts | all exit 0 and all 32 controls fire. Each output equals the committed one line for line, apart from the timings, including every residual in e-notation |
| rerun the certificates | `python3 notes/warp/review_two_scalars_ideals.py --certs --no-controls` reproduces `review_two_scalars_ideals_certs_output.txt`. *Added by the coordinator:* that is the command for that file, which the table above does not state |
| re-solve | #65's D for a generic b(t, x, ρ), fully symbolic, with the coordinator's own Ricci routine (H1) |
| independent check `ind79.py` | run 2 exits 0 in 2 s, and all 9 mutations fire on target. Run 1 stopped at an assertion in the coordinator's code |

`ind79.py` imports nothing from the repository. Its Ricci routine and its polar and Cartesian routes are the
coordinator's, from `I77-mac` and `I78-mac`. This review's cofactors enter it only as data, multiplied into the
coordinator's own coefficients.

| id | what | run 2 |
|---|---|---|
| H0 | planted: Painlevé–Gullstrand Ricci-flat | pass |
| H1 | generic b(t, x, ρ): R(e_ϕ, ·) = 0, the block of §1 and det = D, exactly | pass |
| H2 | c₀ and c₃ as printed, each c_k in (A, B²), and this review's 39 + 39 cofactors give A and B² | pass |
| H3 | 2w = rw′ and w′ + rw″ = 0 force w′ = 0 (`dsolve`) | pass |
| H4 | the accelerating entries at 3 exact points (Cartesian route); the 25 accelerating cofactors give 1 | pass |
| H5 | the ω formula on 6 random instances, uniqueness (rank 4), a violated random dk | pass |
| H6 | b(ρ): R = h((dρ)² − ν′²) and (ν ∧ dν)_{txρ} = b′ | pass |

So D, the ideals by both certificates, the step to "no open set", the accelerating unit ideal, ω and b(ρ) are confirmed
a third time, by code that shares nothing with `two_scalars.py` or with this review's scripts. The lemma's proof (§2)
was read and holds. With this review, #65's status words stand. two_scalars.md now carries the review sentence of §11
and wording points 1, 2 and 4 of §8, as marked sentences.

**Discrepancies.**
- `nonaxi` was silent on its named target `D_cartan_jets`, and the target list was changed after the run. §7 records
  this, and `D_cartan_jets` has five controls that fire.
- Claim A1 was false (`D_slice_nonaxi` 0 of 9). It was replaced by `D_slice_bz0` and `D_slice_offslice` with the
  coordinator's approval. It does not concern #65.
- `acc-drop-odd` was stopped at 12 min, and its method changed to a numeric pre-check; its target did not change.
- Three bugs in this review's own checks failed first runs. They were fixed, with the logs kept.
- The command that writes the certificates file is not stated in the note's table (see the step table above).

**Recommended tasks.** A review of `two_scalars.py`'s algebraic-stage code (`sign_of`'s sampling) is declined. The wall
verdict rests on D ≠ 0 and the lemma, so `split_form` and `sign_of` never decide it. At the sample points, D ≠ 0 was
recomputed here and in §9. The sampled signs of the control's split hold by hand on 0 < β < 1 (§8, item 4).

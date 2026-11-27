# The exchange identities and the Clausen integrals behind Ω₂ and Ω₃, derived again

November 2026. Issue #209 (milestone G5, Paper III), the independent re-derivation asked for by #26's task 4 and #208's E3.
The derivations of §2 to §5 were written from the order equations and the definitions alone, before `omega2.md` and
`omega3.md` were opened; §8 is the comparison made afterwards. The script is
[`exchange_review.py`](exchange_review.py), its output [`exchange_review_output.txt`](exchange_review_output.txt), the raw
logs `I209-rpi5` of negative-mass-raw.

Status words as in AGENTS.md: *derived* (by hand, each step checked), *checked symbolically* (sympy residual 0),
*numerical evidence* (agreement of integrals). A digit count is −log₁₀ of a relative difference between two evaluations.
It is not an error bar.

## Result

TBD

## 1. Setting

Variables on the throat zone. R = r^s = sec²θ, θ ∈ [0, π/2), so that

```math
u = \sqrt{R-1} = \tan\theta,\qquad l = \ln R = -2\ell,\qquad \ell = \ln\cos\theta,\qquad S = \sin^2\theta = 1 - 1/R .
```

The order equations (README.md "Method", `order4.py`'s docstring), for k ≥ 0:

```math
\frac{dY_k}{dl} = F_{k-1} = \sum_{j<k} L_{k-1-j}\,\frac{l^j}{j!},\qquad (R-1)L_k' - \tfrac12 L_k = h_k,\qquad
h_k = (1 - 1/R)L_{k-1} + \sum_{j<k} Y_j c_{k-j} + \tfrac12 Y_k ,
```

with Y₀ = 1, h₀ = Y₀/2, and Y_k = 0 at l = 0 for k ≥ 1 (the mode is normalized by Y = 1 at the throat). The functions c_m are
the table of `cc` (order4.py, order5.py). Used here: c₁ = Ω₀ + 3l/2 − 1 + 1/R, and c_m(l = 0, R = 1) = Ω_{m−1}.
The solution analytic at the throat, and the condition at order k:

```math
L_k = -2\Omega_{k-1} + u\,C_k,\qquad C_k(\theta) = 2\int_0^\theta \frac{h_k - \Omega_{k-1}}{\sin^2\varphi}\,d\varphi,\qquad G_k = C_k(\pi/2) = 0 .
```

**Hypothesis H_k.** G_j = 0 for 2 ≤ j ≤ k − 1 (G₁ = 0 identically). That is, Ω₀ … Ω_{k−3} have their values. Ω_{k−2}
is free, and Ω_{k−1} cancels in G_k. Every identity below is an identity in Ω_{k−2}, and §7 tests it at three values.

**Four facts, used at every boundary.**

- **F1, the substitution.** R − 1 = tan²θ, dR = 2 sec²θ tan θ dθ, so (R − 1)^{−3/2}dR = 2 dθ/sin²θ, which is the second
  form of C_k above. Also dl/dθ = 2 tan θ, hence cot θ · dl/dθ = 2, and tan θ = −dℓ/dθ. *Checked symbolically* (S1, S2).
- **F2, the throat.** h_k(1) = Ω_{k−1}: at R = 1, l = 0, 1 − 1/R = 0, Y_j = δ_{j0}, and c_k(0, 1) = Ω_{k−1}. So
  h_k − Ω_{k−1} = (R − 1)a_k(R) with a_k analytic at R = 1 (induction: l = ln R, 1/R, Y_j, L_{k−1} are analytic there). Then
  C_k = 2∫₀^θ a_k R dφ = O(θ), L_k = −2Ω_{k−1} + O(θ²) is analytic in R − 1 = tan²θ, and Y_k = ∫₀^l F_{k−1}dl′ = O(l) = O(θ²)
  for k ≥ 1. At the throat ℓ = −θ²/2 + O(θ⁴) (S7).
- **F3, growth at θ → π/2** (x = π/2 − θ, cot θ ≈ x, tan θ ≈ 1/x, l ≈ 2 ln(1/x)). Under H_k, every Y_j, L_j, h_j, F_j with
  j ≤ k − 1 (and Y_k) is O(l^p) for some p. Induction: L₀ = −1. If the L_i, i < j, are O(l^p), then Y_i for i ≤ j, F_{j−1}
  and h_j are polynomially bounded in l. Because G_j = 0, C_j(θ) = −2∫_θ^{π/2}(h_j − Ω_{j−1})/sin²φ dφ, and with 1/sin²φ ≤ 2 on
  [π/4, π/2] and l(φ) increasing, |C_j| ≤ K x l^q, so |u C_j| ≤ K l^q. Hence L_j is O(l^q) and Y_{j+1} = ∫F_j dl is O(l^{q+1}).
  Without G_j = 0 the term u G_j ~ G_j/x grows like e^{l/2}, and every integral below diverges. *Numerical evidence*
  (GROW): Y_k ∝ l^{k+1}, L_k ∝ l^k at l = 256.
- **F4, convergence.** By F2 and F3 each integral below converges at both ends, so each split of an integral is legal.
  At θ → 0 the integrands are bounded (O(θ²)/sin²θ), at θ → π/2 they are polynomials in l times at most a bounded function,
  and ∫₀^{π/2} lⁿ dθ < ∞.

## 2. E1: Y_k out of G_k

By F1, (1 − 1/R)/sin²θ = 1, so with S_k = Σ_{j<k} Y_j c_{k−j}:

```math
G_k = 2\int_0^{\pi/2}\frac{h_k - \Omega_{k-1}}{S}\,d\theta
    = 2\int L_{k-1}\,d\theta + 2\int\frac{S_k - \Omega_{k-1}}{S}\,d\theta + \int\frac{Y_k}{S}\,d\theta .
```

Integrate the last term by parts with d(−cot θ) = dθ/sin²θ:

```math
\int_0^{\pi/2}\frac{Y_k}{\sin^2\theta}\,d\theta = \Big[-\cot\theta\,Y_k\Big]_0^{\pi/2} + \int_0^{\pi/2}\cot\theta\,\frac{dY_k}{d\theta}\,d\theta,
\qquad \frac{dY_k}{d\theta} = F_{k-1}\frac{dl}{d\theta} = 2\tan\theta\,F_{k-1} .
```

- **Boundary at θ = 0:** cot θ ≈ 1/θ and Y_k = O(θ²) (F2), so cot θ Y_k = O(θ). Its leading term is F_{k−1}(0)θ with
  F_{k−1}(0) = L_{k−1}(0) = −2Ω_{k−2}.
- **Boundary at θ = π/2:** cot θ ≈ x and Y_k = O(l^p) (F3), so cot θ Y_k = O(x l^p) → 0.
- The integrand is cot θ · 2 tan θ F_{k−1} = 2F_{k−1}.

So ∫Y_k/S = 2∫F_{k−1}dθ, and F_{k−1} = L_{k−1} + Σ_{j=1}^{k−1} L_{k−1−j} l^j/j!:

```math
\boxed{G_k = 4\int L_{k-1}\,d\theta + 2\int\sum_{j=1}^{k-1}L_{k-1-j}\frac{l^j}{j!}\,d\theta + 2\int\frac{S_k - \Omega_{k-1}}{\sin^2\theta}\,d\theta .}
```

This is E1. It holds under H_k, with Ω_{k−2} free. Consequence: Ω_{k−2} enters the right side only through
L_{k−1} = −2Ω_{k−2} + u C_{k−1} (C_{k−1} does not contain it) and through Y₀c_k + Y₁c_{k−1}, where it appears as
Ω_{k−2}(l + Y₁) = 0. So the slope of G_k in Ω_{k−2} is 4·(−2)(π/2) = −4π, the value stated in `order4.py`.

## 3. E2: ∫L_m by a second exchange, m = k − 1

Write L_m = −2Ω_{m−1} + tan θ C_m and use tan θ = −dℓ/dθ, C_m′ = 2(h_m − Ω_{m−1})/sin²θ:

```math
\int_0^{\pi/2}L_m\,d\theta = -\pi\Omega_{m-1} + \Big[-\ell\,C_m\Big]_0^{\pi/2} + 2\int_0^{\pi/2}\frac{\ell\,(h_m - \Omega_{m-1})}{\sin^2\theta}\,d\theta .
```

**Boundary terms of this step.**

- **θ = 0:** ℓ = O(θ²) and C_m = O(θ) (F2), so ℓ C_m = O(θ³).
- **θ = π/2:** this is where G_m = 0 is used. C_m(θ) = G_m − 2∫_θ^{π/2}(…) = −2∫_θ^{π/2}(h_m − Ω_{m−1})/sin²φ dφ = O(x l^q)
  (F3), and ℓ = −l/2, so ℓ C_m = O(x l^{q+1}) → 0. With G_m ≠ 0 the term would be −ℓ G_m, which diverges like ln(1/x), and
  so would ∫L_m (L_m ≈ G_m tan θ).

The remaining integral splits (F4) into 2∫h_m ℓ/sin² − 2Ω_{m−1}∫ℓ/sin². The antiderivative of ℓ/sin²θ is
A(θ) = −cot θ ℓ − θ (d/dθ: ℓ/sin²θ + 1 − 1, since d(−cot θ ℓ) = ℓ/sin²θ + cot θ tan θ). Its boundary values:

- **θ = 0:** cot θ ℓ ≈ −θ/2 → 0, so A(0) = 0.
- **θ = π/2:** cot θ ℓ ≈ x ln x → 0, so A(π/2) = −π/2.

So ∫₀^{π/2} ℓ/sin²θ dθ = −π/2 (S6), and −2Ω_{m−1}(−π/2) = +πΩ_{m−1} cancels −πΩ_{m−1}:

```math
\boxed{\int_0^{\pi/2}L_{m}\,d\theta = 2\int_0^{\pi/2}\frac{h_{m}\,\ln\cos\theta}{\sin^2\theta}\,d\theta .}
```

This is E2 at m = k − 1. It needs H_{m+1} through G_m = 0, nothing about Ω_{m−1}, and it holds for every constant Ω_{m−1} in
h_m. The boundary terms are: [−ℓ C_m] and [−cot θ ℓ], each zero at both ends. #26 checked them numerically only.

## 4. E3: Y_{k−1} out of G_k, with m = k − 1

Apply E1, then E2 at m to 4∫L_m = 8∫h_m ℓ/S. Split h_m = H_m + ½Y_m, with H_m = (1 − 1/R)L_{m−1} + Σ_{j<m}Y_j c_{m−j} free
of Y_m, and S_k = S_k′ + Y_m c₁, with S_k′ = Σ_{j<m}Y_j c_{k−j}. Y_m then enters in two places:

- from E2, through the ½Y_m of h_m: 8∫(½Y_m)ℓ/S = 4∫Y_m ℓ/S;
- from the j = m term of S_k: 2∫Y_m c₁/S.

By F1, c₁ = Ω₀ + 3l/2 − 1 + 1/R = Ω₀ − 3ℓ − S, so c₁/S = (Ω₀ − 3ℓ)/S − 1 (S3). The two terms add to

```math
\int Y_m\Big[\frac{4\ell + 2\Omega_0 - 6\ell}{S} - 2\Big]d\theta = -2\int Y_m\Big[\frac{\ell - \Omega_0}{\sin^2\theta} + 1\Big]d\theta = 2\int Y_m P'\,d\theta,
\qquad P = \cot\theta\,(\ell - \Omega_0) .
```

Here P′ = −[(ℓ − Ω₀)/sin²θ + 1] (S4, S5). Integrate by parts, with dY_m/dθ = 2 tan θ F_{m−1} and
P · 2 tan θ = 2(ℓ − Ω₀) (S4):

```math
2\int_0^{\pi/2} Y_m P'\,d\theta = 2\Big[Y_m P\Big]_0^{\pi/2} - 2\int P\,\frac{dY_m}{d\theta}\,d\theta = 2\Big[Y_m P\Big]_0^{\pi/2} - 4\int_0^{\pi/2}(\ell - \Omega_0)\,F_{m-1}\,d\theta .
```

**The boundary term [Y_m P]** (for k = 5, m = 4, the j = 4 boundary):

- **θ = 0:** ℓ = O(θ²), so P = −Ω₀/θ + O(θ), and Y_m = O(θ²) (F2): Y_m P = O(θ). Its leading term is
  −Ω₀F_{m−1}(0)θ = 2Ω₀Ω_{m−2}θ.
- **θ = π/2:** P = tan x (ℓ − Ω₀) = O(x l), and Y_m = O(l^{m+1}) (F3): Y_m P = O(x l^{m+2}), which is O(x l⁶) → 0 at m = 4.

Both ends vanish, and no Y_m is left:

```math
\boxed{G_k = 8\int\frac{H_m\,\ell}{\sin^2\theta}\,d\theta + 2\int\big(F_m - L_m\big)\,d\theta + 2\int\frac{S_k' - \Omega_{k-1}}{\sin^2\theta}\,d\theta - 4\int(\ell - \Omega_0)\,F_{m-1}\,d\theta,\qquad m = k - 1 .}
```

Also the by-parts identity by itself: ∫Y_m[(ℓ − Ω₀)/sin²θ + 1]dθ = 2∫(ℓ − Ω₀)F_{m−1}dθ. This is E3. It is valid for every
m ≥ 1 (the argument above does not use m = 4), so at k = 4 it removes Y₃ and at k = 5 it removes Y₄.

**At k = 5** (m = 4, F₃ = L₃ + L₂l + ½L₁l² + ⅙L₀l³, F₄ − L₄ = L₃l + ½L₂l² + ⅙L₁l³ + (1/24)L₀l⁴):

```math
G_5 = 8\int\frac{H_4\,\ell}{\sin^2\theta} + 2\int\Big(L_3 l + \tfrac12 L_2 l^2 + \tfrac16 L_1 l^3 + \tfrac1{24}L_0 l^4\Big)
 + 2\int\frac{\sum_{j\le3}Y_jc_{5-j} - \Omega_4}{\sin^2\theta} - 4\int(\ell - \Omega_0)\big(L_3 + L_2 l + \tfrac12 L_1 l^2 + \tfrac16 L_0 l^3\big).
```

Only Y₀ … Y₃ and L₀ … L₃ remain. Y₄, Y₅ and L₄ are gone. The slope of G₅ in Ω₃ is again visible: it comes only from
8∫(Y₀Ω₃)ℓ/S = −4πΩ₃, because Y₀Ω₃l + Y₁Ω₃ = Ω₃(l − l) = 0 in the third term.

```text
G_5 = 0
├─ E1   Y₅ → F₄                       parts on 1/sin²; boundary [−cot θ Y₅]
├─ E2   ∫L₄ → ∫h₄ ℓ/sin²             parts on tan θ; boundaries [−ℓ C₄], [−cot θ ℓ]; needs G₄ = 0
│   └─ h₄ contains ½Y₄
└─ E3   Y₄ → F₃                       parts against P; boundary [Y₄ P]; Y₄ from E2's h₄ and from h₅'s j = 4 term
remains: Y₀ … Y₃, L₀ … L₃
```

## 5. The Clausen integrals

With Cl₂(x) = Σ_{n≥1} sin(nx)/n² (Clausen), C(θ) = Cl₂(π − 2θ) = Σ (−1)ⁿ⁺¹ sin(2nθ)/n², using sin(nπ − 2nθ) =
(−1)ⁿ⁺¹ sin 2nθ. The series converges uniformly, since Σ1/n² < ∞. Cl₂′(x) = −ln(2 sin(x/2)), so
C′(θ) = −2Cl₂′(π − 2θ) = 2 ln(2 cos θ) (S9). C(0) = Cl₂(π) = 0 and C(π/2) = Cl₂(0) = 0.

**∫₀^{π/2} tan θ C dθ = π³/12, route 1 (by parts).** With tan θ = −dℓ/dθ:

```math
\int_0^{\pi/2}\tan\theta\,C\,d\theta = \Big[-\ell\,C\Big]_0^{\pi/2} + \int_0^{\pi/2}\ell\,C'\,d\theta = 2\int_0^{\pi/2}\ell\,(\ell + \ln 2)\,d\theta .
```

Boundary: at 0, ℓ = 0 and C = 0. At π/2, C = Cl₂(2x) = 2x(1 − ln 2x) + O(x³) and ℓ = ln sin x ≈ ln x, so ℓC = O(x ln²x) → 0 (S7).
The moments come from B(a) = ∫₀^{π/2}cos^aθ dθ = (√π/2)Γ((a + 1)/2)/Γ(a/2 + 1) by differentiating in a at a = 0:
∫ℓ = B′(0) = −(π/2) ln 2 and ∫ℓ² = B″(0) = (π/2)(ln²2 + π²/12) (S8). So

```math
2\Big[\tfrac\pi2\big(\ln^2 2 + \tfrac{\pi^2}{12}\big) - \tfrac\pi2\ln^2 2\Big] = \frac{\pi^3}{12} .
```

**Route 2 (the series).** I_n = ∫₀^{π/2}tan θ sin(2nθ)dθ. From sin 2(n + 1)θ − sin 2(n − 1)θ = 2 cos 2nθ sin 2θ and
tan θ sin 2θ = 2 sin²θ = 1 − cos 2θ:

```math
\tan\theta\,[\sin 2(n+1)\theta - \sin 2(n-1)\theta] = 2\cos 2n\theta - \cos 2(n+1)\theta - \cos 2(n-1)\theta .
```

For n ≥ 2 the right side integrates to 0 over [0, π/2], so I_{n+1} = I_{n−1}. Direct integrals give I₁ = ∫2 sin²θ = π/2 and
I₂ = ∫4 sin²θ cos 2θ = −π/2. So I_n = (−1)ⁿ⁺¹π/2 (S10). Termwise integration is legal by Tonelli: with x = π/2 − θ,
|tan θ sin 2nθ|/n² ≤ min(2/n, 1/(x n²)), whose sum over n is O(1 + ln(1/x)), integrable. Then
∫tan C = Σ (−1)ⁿ⁺¹/n² · (−1)ⁿ⁺¹π/2 = (π/2)ζ(2) = π³/12 (S11). The two routes share nothing but the definition of Cl₂.

**∫₀^{π/2} C² dθ = π⁵/360.** The functions sin(2nθ), n ≥ 1, are orthogonal on [0, π/2] with ∫sin²(2nθ) = π/4 (S10).
The series for C converges uniformly, so C² integrates termwise:

```math
\int_0^{\pi/2} C^2\,d\theta = \sum_{n\ge1}\frac{1}{n^4}\cdot\frac{\pi}{4} = \frac\pi4\,\zeta(4) = \frac{\pi^5}{360} .
```

This is Parseval's identity for the orthogonal system sin 2nθ, but only orthogonality and uniform convergence are used, not
completeness.

TBD (checks, third grid, controls, comparison, limits, files)

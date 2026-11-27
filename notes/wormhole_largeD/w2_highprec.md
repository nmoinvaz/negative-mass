# #195: W²(D) of the Ellis wormhole at large D to 60 digits, and Ω₁ and Ω₂ from it

Status: W²a² at 64 values of D from 50 to 10⁴, and at growth.csv's 16 values, is **numerical evidence** to 61 digits (two lanes
of one method agree). Ω₁ = 2π²/3 + 8 ln²2 is confirmed by the fit to 8e-50 (the median; the worst cell, 1.7e-45). Ω₂ = 74.36938672793331759538467029584661726673414
is **identified** as 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 by PSLQ (an identification, not a derivation; #26 asks for the derivation).
Every number below is in `w2_output.txt` (the run of record; raw folder `negative-mass-raw/I195-dh610`, README first, last commit
when the output was assembled in its first line) unless it names another raw file. No sentence here claims novelty.

## Result

- **W²(D).** `w2_crate` solves E2 of `master.py` in f64, double-double, quad-double and octo-double. The quad-double and
  octo-double values agree to 60.8 to 64.6 digits (median 62.4) at all 64 values of D of the list, 50 to 10⁴. Two lanes of the
  same method agree; that is an estimator of the error, not a bound (what the lanes share is under "What the checks are").
- **growth.csv.** The eight `modes_r.py` rows (D ≥ 20) are low by 8e-10 at D = 20, growing to 1.9e-6 at D = 1000 and 6.0e-6 at
  D = 2000. The criterion G1 (every row within half a unit of its 12th decimal) is **not met** in 16 of 16 rows. The D = 4 … 11
  rows (TS13) differ from the crate by 9e-12 to 6.2e-10. `modes_r.py`'s own tighter run (X 50, ρ₀/10, rtol 1e-13, in
  `logs_modes_r.txt`) moves each row toward the crate's value, by 42 % (D = 20) to 84 % (D = 2000) of the gap.
- **Ω₁** = 10.42336037873851714322648087719742 (F3k's median over 6 admitted fit cells, range 1.7e-45). Against
  2π²/3 + 8 ln²2 the median differs by 8.0e-50; the six cells differ from it by 1.1e-50 to 6.7e-48, and one by 1.7e-45. PSLQ finds 3Ω₁ − 2π² − 24 ln²2 = 0 at residual 2.4e-49, coefficients ≤ 100, and fails for a perturbed basis. The
  earlier identification was at 1e-9.
- **Ω₂** = 74.369386727933317595384670295846617267 (median, range 3.2e-40: 41 digits as an estimator range). PSLQ in
  {Ω₂, ζ(3), π² ln 2, ln³2} finds (3, −240, 8, 32) at residual 2.7e-44 (45 digits of Ω₂), and fails for a perturbed basis.
  The identified closed form 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 agrees with the fit's median to all 45 digits printed
  (`dev_27_omega2_vs_closed.log`). The README's "about 73" was a numerical estimate; Ω₂ = 74.37.
- **Ω₃** = −652.8401075700561254 (range 2.7e-35) and **Ω₄** = −12794.49986002 (range 1.5e-30) from the same cells. A try at Ω₃ in the
  pure weight-4 constants found no relation (`dev_25_pslq_omega3.*`, a development run, not a result).
- **Not met, and why** (details below): S4, S4w, S5 (planted checks), G1, R4, F3. S5s, E0, M1, O2, L1, F1, F2, F3k, P1, P2, P3
  pass. The driver exits 1 by design (14 of 20 checks pass).

In the repository's form, with s = 2(D − 3) and ε = 1/s, the series is

```math
W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\pi^2/3 + 4\ln^2 2}{D-3}
 + \frac{20\zeta(3) - \tfrac23\pi^2\ln 2 - \tfrac83\ln^3 2}{(D-3)^2} - \frac{81.60501344625702}{(D-3)^3} + O((D-3)^{-4}),
```

the coefficients 5.21168018936925857, 18.59234668198332940 and (fitted) −81.605013446257016. The first two are the identified
closed forms; the last is the fit.

## What it is

- **Crate** `w2_crate/` (Cargo.toml, Cargo.lock, src/); the driver `w2_highprec.py` builds nothing: `cargo build --release -j 2`
  in the crate first. `python w2_highprec.py --logs=<crate logs>` runs the checks and the fits. `--mutate=NAME` runs a control.
- **The equation.** ℓ = s ln r, r = e^{εℓ}, F = 1 − e^{−ℓ}, w = dY/d ln r, W² = s + Ω:

  Y_ℓ = ε w,  F w_ℓ = A Y + B w,  A = (1 + εΩ)e^{2εℓ} − (½ + ε) + ε e^{−ℓ},  B = (½ + 2ε) − 2ε e^{−ℓ}.

  It is E2 of `master.py` (check E0, symbolic: residual 0) and `modes_r.py`'s and `matched.py`'s system in other variables.
- **Throat.** ℓ = 0 is a regular singular point with exponents 0 and ½ (0 and 3 in proper distance). The even mode Y(0) = 1 is
  the analytic one: w₀ = −(1 + 2εΩ), a recursion with denominators (m − ½), radius of convergence 2π (zeros of F at 2πik);
  summed at ℓ₀ = 1. Its first three orders equal `matched.py`'s Y₀, Y₁, Y₂ (check O2: the difference at D = 10⁴ is
  −9.8 ε³, an O(ε³) term).
- **Far end.** The two solutions grow like r and rⁿ, the gap of their exponents in ℓ is at least ½. The linear system is
  integrated inward from ℓ_∞ (190 to 325 in quad-double, 255 to 623 in octo-double) from an arbitrary vector; the other
  solution falls by e^{−T}, T = 162 (qd), 311 (od). No Bessel function and no truncation of F.
- **Eigenvalue.** Newton on M(Ω) = Y_in w_ser − w_in Y_ser at ℓ₀, the exact derivative carried as dual parts of the recursion;
  5 or 6 iterations in qd and 6 to 8 in od, from the three-term series (or growth.csv's value for D < 50). Taylor stepping by multishoot's `plan_path`
  rule (the largest h with the last two terms below tol × the state), in τ = (ℓ − ℓ₀)/h so coefficients stay near 1.
- **Lanes** (f64, dd, qd, od behind one trait, as #186's crate; reciprocals 1/k formed in the lane, #193's lesson; the control
  `recip_f64` rounds them from f64 and fails S1, S2, S3, S5s and F3k):

  | lane | Taylor order | tol | series order | T | stated digits (planted) |
  |---|---|---|---|---|---|
  | f64 | 24 | 1e-17 | 20 | 51 | 11 |
  | dd | 48 | 1e-33 | 46 | 88 | 27 |
  | qd | 96 | 1e-65 | 90 | 162 | 56 |
  | od | 180 | 1e-130 | 170 | 311 | 110 |

- **Libraries.** The crate depends on multishoot 4.2.0, multisolver 0.2.2 and multidouble (the pin of `taylor_polar`). *Where
  they fit:* multishoot's `Equations` has what this problem is (a Frobenius endpoint in `left_state`, an eigenvalue as a global,
  inward segments in `Shooting`), and `Plan::optimal_kt` is used as a cross-check of the order (45, 83, 159 for dd, qd, od
  against 48, 96, 180 chosen here, printed as `kt_ms`). *Where they did not:* `Equations` takes its constants in
  double-double (#193), so a quad- or octo-double lane would be double-double accurate in its Taylor recursion; it has no f64 lane;
  and the problem has one unknown, so multisolver's `newton` (a `Problem` over part arrays) would add a layer to a ten-line
  Newton step with an exact derivative. The step rule of `plan_path` is therefore re-implemented, the engine is not called, and
  multisolver is not called. No change to a library was needed or made; there is no pull request. A multiple-shooting second
  method (`Shooting` with nodes) was not built; the digits are checked by the lane and by the independent comparisons below.

## What the checks are, and what they share

Defined in the raw README (A1, A2 for what was added by rulings). The lane check (L1) compares qd with od: the same
equation, series, flow and Newton code, so it catches arithmetic and truncation (T differs: 162 against 311) but not a wrong
equation. Independent of that: E0 (sympy, from E2's text), M1 (`master.py`'s scipy integration in proper distance, D = 5 … 8,
worst difference 6.1e-4 against its stated 1e-3), TS13's printed values (D = 20, 50, 100 agree with the crate to 1.4e-10,
7e-11, 6.6e-10, those digits as README.md quotes them), and the planted eigenproblem.

### Planted eigenproblem (same code; A is replaced)

A_pl = ε(Ω − Ω*) + sκ u [F(1 + 2κu) + B], u = e^{−ℓ}/(1 + κF), κ = min(ε, 10⁻³), Ω* = 5/2. At Ω* its solution is
Y* = 1/(1 + κF), so the eigenvalue is exactly Ω*. D = 5, 50, 1000, 10⁴; log10 of the relative error (stated digits above):

| lane | D | S1 series | S2 flow | S3 contraction | S4 Ω, start +0.3 | S4 Ω, start −0.4 | S4w (W²), start −0.4 | S5 | S5 bound |
|---|---|---|---|---|---|---|---|---|---|
| f64 | 1000 | −14.8 | −15.3 | −15.8 | −10.56 | −10.37 | −13.27 | −10.52 | −11.0 |
| f64 | 10⁴ | −16.0 | −15.0 | −15.8 | −9.02 | −9.94 | −13.84 | −8.69 | −11.0 |
| dd | 1000 | −31.0 | −31.6 | −31.9 | −26.84 | −27.53 | −30.43 | −23.14 | −25.77 |
| dd | 10⁴ | −31.4 | −31.3 | −32.1 | −25.20 | −26.25 | −30.15 | −21.09 | −25.77 |
| qd | 10⁴ | −64.7 | −63.7 | −64.1 | −59.11 | −57.47 | −61.37 | −37.16 | −41.86 |
| od | 10⁴ | −129.5 | −129.1 | −129.9 | −123.74 | −123.49 | −127.39 | −69.52 | −74.22 |

All 16 (lane, D) cases are in `w2_output.txt`. **S1, S2, S3 (series, flow, contraction): met in 16 of 16.** **S4 as worded
(Ω, both starts): met in 21 of 32, not met** in f64 and dd at D = 1000 and 10⁴ (rows above), and for the start Ω* − 0.4 at
D = 5 in every lane. **S4w (added by ruling: W² relative, the same digits): met in 28 of 32**, failing only at that D = 5
start. Ω enters A only through εΩ, so Ω's condition number is about 1/ε = s; the fit uses qd, where S4 as worded is met at
D ≥ 50. **S5 as worded (the shorter contraction, bound 100e^{−T′}): met in 8 of 16, not met** in every lane at D = 1000 and
10⁴. **S5s (100 s² e^{−T′}): met in 12 of 12** (dd, qd, od), s from Ω's conditioning and s from the other solution's w/Y.

**The start below Ω\*.** From Ω* − 0.4 at D = 5, Newton reaches no root in dd and qd (60 iterations, Ω = −1.26 and −0.21,
drifting 0.09 per step on a flat positive M), and f64 stops at Ω = −1.5041. Below Ω ≈ −1.5, M(Ω) oscillates with a zero every
0.025 to 0.05 and the node count of Y rises from 0 to 22 at Ω = −4: those are the planted system's higher eigenvalues
(eigenfunctions with nodes), and f64 lands on the first (raw `dev_10` to `dev_13`). The same start at D = 6 also fails (Ω = −3.509 in dd and f64, `dev_10`, `dev_11`); at D = 8 and 11 it reaches Ω* (f64, `dev_11`).

### The real problem's start rule (R4, not met)

Newton from the start shifted by +0.5 reaches the same root to 1e-50 at D = 4 … 11, 50, 1000. Shifted by −0.5 it does **not
converge** at D = 4 … 11 (40 iterations, W² between 0.24 and 2.6 and still moving) and reaches no other root; at D = 50 and
1000 both shifts agree. So the start rule must start above the root or close to it, as the record's runs do (growth.csv's
value for D < 50, the three-term series otherwise, both within 1e-4 of the root).

## growth.csv, and the tighter `modes_r.py` run

G1 is not met in 16 of 16 rows (|crate − csv| ≥ 9e-12). The comparator passes the crate's own values printed to 12 decimals, so
the failures are the csv's. `logs_modes_r.txt` was not rerun.

| D | source | csv | crate (qd) | crate − csv | tighter − csv | fraction of the gap closed |
|---|---|---|---|---|---|---|
| 4 | TS13 | 1.397052433715 | 1.39705243433350684 | 6.19e-10 | | |
| 10 | TS13 | 12.044265014744 | 12.044265014497883 | −2.46e-10 | | |
| 20 | modes_r | 31.575110127562 | 31.5751101283584884 | 7.96e-10 | 3.38e-10 | 0.42 |
| 50 | modes_r | 91.345775903365 | 91.345775913629638 | 1.03e-8 | 7.13e-9 | 0.70 |
| 100 | modes_r | 191.283017684226 | 191.28301773035896 | 4.61e-8 | 3.57e-8 | 0.77 |
| 200 | modes_r | 391.254334217002 | 391.25433438080667 | 1.64e-7 | 1.30e-7 | 0.79 |
| 500 | modes_r | 991.237971426592 | 991.23797214792929 | 7.21e-7 | 6.05e-7 | 0.84 |
| 1000 | modes_r | 1991.232655318753 | 1991.2326572612893 | 1.94e-6 | 1.46e-6 | 0.75 |
| 2000 | modes_r | 3991.230019738767 | 3991.2300256842525 | 5.95e-6 | 5.02e-6 | 0.84 |

The tighter run moves toward the crate at all eight rows and does not reach it. The rows D = 5 … 9, 11, 30 are in
`w2_output.txt`. The rows D = 4 … 11 are labelled TS13 Table I in the csv; README.md quotes TS13's table to 8 to 10 digits, and
the csv's 12 decimals differ from the crate by up to 6e-10, so the csv's last digits do not all come from TS13's table.
*Update from #202:* that inference is wrong. TS13's Table I prints 15 digits (arXiv:1309.2058, the TeX source), and the csv's
D = 4 … 11 rows are those values rounded to 12 decimals. The 9.9e-12 to 6.4e-10 gap to this crate is TS13's own error:
`one_mode.log`'s independent W(r) agrees with the crate to 12 digits (`ts13_digits.py`, notes/literature.md).

## W²(D) at the list

The 64 values, 60 digits each, with the digits shared by qd and od, are in `w2_output.txt`. A selection (qd):

| D | W²a² | digits shared with od |
|---|---|---|
| 50 | 91.34577591362963760699372803881622901592398414253952989481 | 63.9 |
| 98 | 187.28422641050103762611373955609994282896051671148215751305 | 63.8 |
| 192 | 375.25549408559253011103199827743964968933823427491764366 | 63.3 |
| 484 | 959.23832598428817942794497348203608338336373203053518573566 | 62.4 |
| 1032 | 2055.2324935622603317264158516121507570762480475695439807121 | 62.5 |
| 2023 | 4037.2299958640103093865819044685247359076253747315392778160 | 61.5 |
| 5103 | 10197.228433890035526247897268186293859846914096223681230862 | 61.3 |
| 10000 | 19991.227932788129793135221624472658961012666469710153134137 | 61.2 |

The loss from 64 digits (quad-double) to 61 is the conditioning: W² = s + Ω with Ω's absolute error about s × the lane's
roundoff. The qd–od difference over the 64 values has median 62.4 shared digits and minimum 60.8.

## The fit

Data: the 64 qd values, g = W²a² − 2(D − 3) + 4 ln 2 = Ω₁ε + Ω₂ε² + …; least squares of g/ε on ε^{k−1}, k = 1 … K, by QR in
mpmath at 300 digits, over cells (K, D ≥ D_lo). The gate admits a cell when its largest relative residual is at most 100 times the
lane disagreement relative to g (δ = 2.97e-54, gate 2.97e-52). The headline is the median over admitted cells and the error bar the
range over them: a spread over fit choices, an estimator range, not a bound.

- **Planted first (F1, F2; met).** Synthetic g from Ω_k = A (k+1)! c^k with Ω₁ = 10.4233603787385, Ω₂ = 73, relative noise 1e-40:
  4 admitted cells, Ω₁ off by 3.2e-37 (range 8.4e-36), Ω₂ by 2.6e-32 (range 1.2e-30). With 0.5 ε³ ln ε added, M1 admits none of 28
  cells and M2 (with log terms) admits one and recovers Ω₂ to 3.3e-24.
- **F3 as declared (K ≤ 24): not met.** One admitted cell (D_lo = 500, K = 24, residual 2.9e-54). A single cell has range 0, so
  the clause "|Ω₁ − closed form| ≤ 3 × range" cannot hold: **a defect of the criterion's wording, not of the data**; the
  difference from 2π²/3 + 8 ln²2 is 5.1e-49 beside it.
- **F3k (met).** M1's grid with K = 28 and 32 added, **after a development run** (dev_21), when the residuals were seen to plateau
  at K ≥ 20; the gate, windows, K + 4 rule, tolerances and headline rule are unchanged. 6 cells are admitted: (D_lo, K) =
  (100, 32), (200, 28), (200, 32), (500, 24), (500, 28), (500, 32). Ω₁ = 10.4233603787385171432265 and Ω₂ = 74.3693867279333175953847
  to the digits shown in every one of them; the ranges are 1.7e-45 and 3.2e-40, five cells agree with the closed form for Ω₂ to 1e-42 and one deviates by 3.2e-40 (that deviation is the range). The headline Ω₁, Ω₂, Ω₃, Ω₄ above are F3k's. Residuals of the other cells fall with K down to 1e-54, the
  data's floor: the series model (M1, no logarithms) describes the data, and M2 was not needed.

## PSLQ (identifications; mpmath, as declared)

| basis | tolerance, H | result |
|---|---|---|
| {Ω₁, π², ln²2} | 1e-11, 100 | found (−3, 2, 24), residual 2.4e-49; budget 7.9 digits against 12.0; the perturbed basis (ln²2 + 10⁻³) finds none |
| {Ω₂, ζ(3), π² ln 2, ln³2} | 1e-13, 1000 | found (3, −240, 8, 32), residual 2.7e-44; budget 14.2 against 14.9 at the declared τ; the perturbed basis (ζ(3) + 10⁻³) finds none |
| {Ω₂, ζ(3), π² ln 2, ln³2, π², ln²2, ln 2, 1} | 3.2e-39, 100 | not found at H = 100 (the relation above has the coefficient 240); budget 19.4 against 40.4 |

The budget is the digits a chance relation of that size would need: for P2, 4 log10(2·1000 + 1) + 1 = 14.2, against the 45
digits to which the relation holds (residual 2.7e-44 on Ω₂ = 74): the relation exceeds the budget by about 31 digits. These
are identifications. With ζ(3), π² ln 2 and ln³2 only, and no lower weights, Ω₂ is of pure weight 3 as Ω₀ is of weight 1 and
Ω₁ of weight 2; the derivation is #26's.

## Controls

Each check has a control; the targets were computed in development runs (raw README A2). A control fires when every target
fails with it. The checks that already fail without a control (S4, S4w, S5, F3) are not targets.

- `wrong_B`: E0, S1, S2, S3, S5s, O2, M1, F3k. `wrong_A_sign`: E0, O2, M1, F3k (the planted A replaces the real one, so S1 to S3
  cannot see it). `recip_f64`: S1, S2, S3, S5s, F3k. `short_series`: S1, S5s, F3k. `step_loose`: S2, S3, S5s, F3k.
  `wrong_planted`: S1, S2, S3. `wrong_eps`: F3k only (the planted Ω term vanishes at Ω*).
- `short_contraction` (T/4): S3, F3k, and S5s in **7 of its 8 cells**; it is **silent at dd D = 5** (−27.33 against −27.0, the
  bound's). T/8 would fire there (−22.56); the control was kept as declared.
- `fit_wrong_eps`: F1, F2, F3k, P1. `fit_drop_term`: F1, F2; **silent on F3k** (3 of 6 cells still admitted, because K = 28
  and 32 leave enough terms). `pslq_wrong_basis`: P1; P2 flips from found to not found.
- The raw folder holds the runs (`dev_19`, `dev_23`). `csv_shift` was dropped: G1 fails on 16 of 16 rows unmodified.

## Timings (dh610-fedora, 20 cores shared with choptuik-large-d, load average 9 to 19 from others, one thread per process, at most 4 at once)

| run | wall | cpu inside the crate | peak RSS |
|---|---|---|---|
| qd, 64 values of D | 50 s | 50.4 s (2106 steps, 1087 retries); 0.6 to 1.4 s per D | 3.5 MiB |
| od, 64 values of D (3 processes) | 16.9, 18.6, 23.3 min | 3526 s in all (2717 steps, 611 retries); 30 to 202 s per D | 3.3 to 3.7 MiB |
| planted, qd (D = 5, 50, 1000, 10⁴) | 1.3 min | D = 5 takes 64 s (a start that does not converge) | 3.5 MiB |
| planted, od (2 processes) | 20.6 and 27.2 min | 276 to 996 s per D | 3.4 to 3.7 MiB |
| the driver | 74 s | | 117 MiB |

The same qd list took 75 s at the development stage (`dev_16`). Memory is reported only.

## Defects and discrepancies found

- **The issue's series** was written in 1/(D − 3) while asking for Ω₁ = 2π²/3 + 8 ln²2; the repository's Ω_k are the coefficients of
  ε^k, ε = 1/(2(D − 3)). The coordinator confirmed and corrected the issue.
- **growth.csv**: all 16 rows differ from the crate by more than the half unit of the 12th decimal (9e-12 to 6.2e-10 for the eight
  TS13 rows, 8e-10 to 6e-6 for the eight `modes_r.py` rows, which are low). README.md's table for D = 20 … 1000 quotes the same values.
- **The first planted construction** (Y* = exp(−κ(1 − e^{−ℓ})), κ = 3/4, Ω weighted by e^{2εℓ}) was a bad test, not a bad
  method: its Taylor coefficients fall by only 0.35 per order, and the real system's contraction rule does not apply to it. The
  coordinator accepted its replacement (A1); next time the replacement is asked about first.
- **Criteria that fail as worded** (S4, S5, F3, and R4 for the start below the root): reported, not changed.
- **A control silent on its target**: `fit_drop_term` on F3k; `short_contraction` on S5s at dd D = 5; `wrong_eps` on every planted check.

## What README.md's "Result" and "Limits", #25 and #26 should now say

- **README.md, Result.** The third term is **computed** (not only to ten digits): Ω₁ = 10.42336037873851714322648088 from a fit of
  64 quad-double values of W², equal to 2π²/3 + 8 ln²2 to 1e-45, and **identified** by PSLQ at 2.4e-49. The fourth term is added:
  Ω₂ = 74.3693867279333175953847 (not "about 73"), **identified** as 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 at 2.7e-44. The table's column
  "W²a² here" for D ≥ 20 should be replaced by this note's values (the printed ones are low by 8e-10 to 6e-6). The sentence that
  the r-integration's values move by 4e-10 to 1.5e-6 under a change of parameters understates their error.
- **README.md, Limits.** "The third term's closed form is an identification, not a derivation" stays true and now holds for the
  fourth term too; what changed is the precision (1e-9 to 1e-49 for Ω₁).
- **#25** (derive Ω₁): this gives an independent number, Ω₁ − (2π²/3 + 8 ln²2) = 8e-50 (F3k, median of 6 cells); it does not replace the
  derivation.
- **#26** (Ω₂): the identified target is 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 = 74.3693867279333175953846703, to be derived from the order-4 integral. Its
  acceptance "Ω₂ to 20 digits": the fit gives 40 (range) and the closed form agrees to 45. "The PSLQ relation fails for a
  perturbed basis": it does. **"The four-term series against growth.csv at D = 200 … 1000 leaves a remainder falling like D⁻³"
  cannot be tested against growth.csv:** its error at D = 1000 and 2000 (1.9e-6, 6.0e-6) is a tenth of the fourth term at D = 1000
  (Ω₂/s² = 1.9e-5) and larger than it at D = 2000 (4.7e-6), and 24 and 580 times the remainder (Ω₃/s³ = −8.2e-8, −1.0e-8). Against the crate's values the remainder
  times s³ is −685.0, −665.7, −659.2, −656.0 at D = 200, 500, 1000, 2000, tending to Ω₃ = −652.84; against growth.csv it is −695, −1374,
  −16060, −379457.

## Recommended tasks

1. **Regenerate growth.csv from the crate.** What: replace the `modes_r.py` rows (D ≥ 20) by this crate's values to 12 decimals, and
   recheck the numbers of Paper III and `profiles.py`, `paper_numbers.py` that read it. Why: the eight `modes_r.py` rows are off by 8e-10 to 6e-6 and the TS13 rows by up to 6e-10.
   Accept: each printed row equals the crate's value to its digits, a check fails on a row changed in its last digit, and the
   papers' quoted growth rates are rechecked. Depends: this issue.
2. **Ω₃ in weight-4 constants.** What: Ω₃ = −652.84010757005612543524 (35 digits from the same fit); get 60 digits (a denser D list, od
   data, K up to 40) and search {Ω₃, π⁴, π² ln²2, ln⁴2, ζ(3) ln 2, Li₄(½)} and the lower weights. Why: Ω₀, Ω₁, Ω₂ are pure weight 1, 2, 3; the
   pure weight-4 basis found nothing at 35 digits. Accept: the relation holds at 1e-50, fails for a perturbed basis, and the budget is met, or
   "not found" is reported with the basis. Depends: this issue.
3. **`Equations` for this problem in multishoot.** What: write the system as a multishoot `Equations` (Frobenius `left_state`,
   Ω as the global, inward segments) once #193's lane-precision constants land, and compare with this crate. Why: a
   second solver of the same equation (L1 shares its code). Accept: qd W² at D = 50, 1000 equal to this crate's to 1e-55, and a
   mutated equation fails. Depends: #193.
4. **TS13's printed digits.** What: read TS13 Table I and compare its printed digits with growth.csv's D = 4 … 11 rows. Why: the csv's
   12 decimals differ from the crate by up to 6e-10 while README.md quotes 8 to 10 digits. Accept: each row's source digits recorded in
   `notes/literature.md`. Depends: none.

## Files

- `w2_crate/` the crate; `w2_highprec.py` the driver; `w2_output.txt` its output of record; this note.
- Raw: `negative-mass-raw/I195-dh610` (README with the criteria and amendments A1, A2; `dev_*` the development runs; `run_*` the runs of
  record, `run_08_driver_final` the output above).

## Coordinator's review

Under review: `wip/i195-w2-largeD` at 514ac42, raw `I195-dh610` at e70f58a. Review folder
`negative-mass-raw/I195-mac-review`, whose README holds the criteria, committed before the first run of `ind195.py`. The
changed files are those the brief names.

- **Rerun on the Mac** of every run of record, with paths mapped from the `.env` files.
  - The 68 od values are bit-identical to dh610's.
  - 4 of the 64 qd values of the D list differ, by 3e-64 to 3e-63, and a few planted fields and Newton traces differ.
    *Added by the coordinator:* the cause is the step planner's f64 `powf` (`flow.rs:161`) and the ℓ_∞ rule's f64 `exp`
    (`flow.rs:191`). Their results differ by an ulp between macOS's libm and glibc, so step sizes differ. multishoot's
    `kth_root` would make the planner bit-portable. Every difference is at qd roundoff, far below every bar.
  - The driver rerun gives the same verdicts, fit and PSLQ numbers.
- **ind195**, run 1, 4 of 4:
  - K1: the coordinator's own transformation of E2 into the ℓ system gives the A and B used here (sympy residual 0).
  - K2: an independent W², in Python and mpmath at 50 digits. It matches at ℓ_m = 3 instead of the throat series at
    ℓ₀ = 1, and starts the far side on the frozen decaying eigenvector at ℓ_∞ = 180 instead of (1, 1). It agrees with the
    crate to 5.6e-51 at D = 50 and 4.2e-40 at D = 1000.
  - K3: Neville extrapolation instead of least squares. It gives Ω₁ within 4e-43 of 2π²/3 + 8 ln²2, and Ω₂ within 1e-39
    of 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 at N = 16, 20 and 24.
  - K4: PSLQ returns (3, −240, 8, 32), and nothing for a perturbed basis.
  - The mutations `k1-B`, `k2-eps` and `k3-omega1` fail exactly their targets.
- **Not met, as ruled:** S4, S5, F3 (one admitted cell), G1 (0 of 16 csv rows) and R4. R4 is new: Newton from the start
  shifted by −0.5 does not converge at D = 4 … 11, and from +0.5 it reaches the same root. These stand as the runner
  reports them. Silent controls (`fit_drop_term` on F3k, `short_contraction` at dd D = 5, `wrong_eps` on the planted
  checks) are reported as silent.
- *Added by the coordinator:* G1's failure is a finding about growth.csv, whose `modes_r.py` rows are low. The
  `logs_modes_r.txt` run at tighter settings moves each row toward this crate's value. The README's table and its
  "66 ± 6" rested on those rows. README.md's "Result" and "Limits" carry an update from this issue.
- Follow-ups: the recommended tasks 1, 2 and 4 are filed. Task 3 waits for a multidouble release with the correctly
  rounded reciprocals (multidouble PR #188; #193 was closed for it).

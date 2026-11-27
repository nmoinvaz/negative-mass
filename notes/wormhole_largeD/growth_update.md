# #200: growth.csv regenerated from #195's crate

Status: **checked** by `check_growth.py` (G1 to G6, R0 to R2 pass; 56 mutations, all on their declared targets). The
D ≥ 20 rows equal #195's quad-double values rounded to 12 decimals. The D = 4 … 11 rows are byte-identical to 613e866
and still differ from the crate by 9.5e-12 to 6.2e-10 (#202). The coordinator applies the text below; nothing outside
the files this issue names was edited. Branch `wip/i200-growth-csv`, raw `negative-mass-raw/I200-rpi5` (README first).

## Result

| acceptance | verdict |
|---|---|
| Each printed row equals the crate's value to its digits | **met for the 8 rows D ≥ 20** (G2, 8 of 8). **Not met for the 8 TS13 rows** (kept as the task states): they equal TS13's printed values (G4, 8 of 8) and differ from the crate (`logs_check_growth.txt:11-18`). |
| A check fails on a row changed in its last digit | met: `last_digit:D` fires G2 (D ≥ 20) or G4, G6 (D ≤ 11) for all 16 rows |
| Each quoted number is listed with old and new value | below; changed where its printed digits change |

Old "73.0" (D = 500) becomes **73.7**, old "66 ± 6" (D = 1000) becomes **74.0** (74.0388 ± 2e-6). Both rested on the
low `modes_r.py` rows: those were low by 7.2e-7 and 1.9e-6, times s² = 988036 and 3976036, which is 0.71 and 7.7. The
old interval 60 to 72 did not contain the new value.

## What was done

```
notes/wormhole_largeD/
├── profiles.py            D >= 20 rows: crate_w2() reads w2_output.txt's table, round12() half-even on the printed digits, source
│                          "w2_crate qd (#195)"; --growth-only skips the profile integrations; the D <= 11 rows are written as before
├── growth.csv             lines 10 to 17 changed (the D >= 20 rows); lines 1 to 9 unchanged
├── check_growth.py        G1-G6 rows, R0-R2 remainder, reported tables, --mutate=NAME (new)
├── logs_check_growth.txt  its output (new; not in the issue's file list, remove if unwanted)
└── logs_paper_numbers.txt paper_numbers.py rerun, 0 checks failed
```

The source string is `w2_crate qd (#195)`, at the coordinator's request: the task's `w2_crate qd, #195` made
`w2_highprec.py`'s `read_csv` (a plain `split(",")`) raise `ValueError`. Tested; `read_csv` reads the new file (16
rows).

## growth.csv, old and new (source: `logs_check_growth.txt:68-75`, `growth.csv:10-17`)

| D | old csv (`modes_r.py`) | new csv (`w2_crate qd (#195)`) | new − old |
|---|---|---|---|
| 20 | 31.575110127562 | 31.575110128358 | +7.96e-10 |
| 30 | 51.440526078330 | 51.440526081063 | +2.73e-09 |
| 50 | 91.345775903365 | 91.345775913630 | +1.03e-08 |
| 100 | 191.283017684226 | 191.283017730359 | +4.61e-08 |
| 200 | 391.254334217002 | 391.254334380807 | +1.64e-07 |
| 500 | 991.237971426592 | 991.237972147929 | +7.21e-07 |
| 1000 | 1991.232655318753 | 1991.232657261289 | +1.94e-06 |
| 2000 | 3991.230019738767 | 3991.230025684252 | +5.95e-06 |

The old rows did not reproduce across machines: the unmodified `profiles.py` on rpi5 gives rows that differ from the
committed ones by 1e-12 (D = 50) to 2.7e-8 (D = 2000), and `profiles.csv` by up to 1.8e-9 relative (raw `dev_01_*`).
The new rows do not depend on the machine. The closest a row comes to a rounding tie is 0.012 of a unit of the 12th
decimal (D = 20), far above the crate's error.

**Independent support for the crate.** `one_mode.log:29,31,33,35` (an unrelated shooting code that shares only the
master equation) gives −E = 1.39705243433, 12.0442650145, 31.5751101284 and 91.3457759136 at D = 4, 10, 20, 50. These
equal the crate's values at every printed digit (`logs_check_growth.txt:11,17,19,21`). At D = 4 and 10 they do not
equal TS13's printed digits (1.39705243372, 12.0442650147). That is evidence for #202.

## paper_numbers.py, rerun on the new csv

`logs_paper_numbers.txt`, old at 613e866. Exit 0, "checks failed: 0", as before. Every number the rerun prints that
changed:

| line | quantity | old | new |
|---|---|---|---|
| 14 | D = 20, W²a² | 31.5751101276 | 31.5751101284 |
| 15 | D = 30, W²a² | 51.4405260783 | 51.4405260811 |
| 16 | D = 50, W²a² | 91.3457759034 | 91.3457759136 |
| 17 | D = 100, W²a² | 191.2830176842 | 191.2830177304 |
| 18 | D = 200, W²a² | 391.2543342170 | 391.2543343808 |
| 19 | D = 500, W²a² | 991.2379714266 | 991.2379721479 |
| 20 | D = 1000, W²a² | 1991.2326553188 | 1991.2326572613 |
| 21 | D = 2000, W²a² | 3991.2300197388 | 3991.2300256843 |
| 19, 20, 21 | relative difference from 3 terms, D = 500, 1000, 2000 | 7.45e-08, 8.38e-09, 3.24e-10 | 7.53e-08, 9.35e-09, 1.17e-09 |
| 18, 19, 20, 21 | (data − 3 terms) s², D = 200, 500, 1000, 2000 | 72.61, 72.99, 66.32, −20.64 | 72.63, 73.70, 74.04, 74.21 |
| 14 to 21 | source column | `modes_r.py` | `w2_crate qd (#195)` |

Unchanged at the printed digits: the D = 20, 30, 50, 100 relative differences (1.30e-03, 3.91e-04, 8.19e-05,
9.82e-06), the (data − 3 terms) s² at D = 20, 30, 50, 100 (47.55, 58.58, 66.07, 70.67), D = 4 … 11, the e-folding
table (D = 100: 0.0723), the ST17 slopes and sigmas, the SH2002 comparison. `paper_numbers.py`'s own checks: none
fails. Stale text in the script, not edited (the task allows edits only for a reader that needs the new source
column):

```diff
 Inputs: growth.csv (profiles.py: Torii & Shinkai's Table I for D <= 11, modes_r.py for D >= 20); logs_order3.txt
+Inputs: growth.csv (profiles.py: Torii & Shinkai's Table I for D <= 11, #195's quad-double values from w2_output.txt for D >= 20); logs_order3.txt
-print("  (the r-integration's own uncertainty, logs_modes_r.txt, times s^2: 0.02 at D = 200, 0.6 at 500, 6 at 1000, 80 at 2000; so the last column means nothing beyond D = 1000)")
+print("  (D >= 20: the csv's rounding to 12 decimals times s^2 is at most 8e-8 at D = 200 and 8e-6 at D = 2000; the column tends to Omega_2 + Omega_3/s)")
-check(all(70 < v < 76 for v in fit.values()), "remainder ~ 73 / s^2 at D = 200, 500")
+check(all(70 < v < 76 for v in fit.values()), "remainder x s^2 between 70 and 76 at D = 200, 500")
```

The removed print line is printed in `logs_paper_numbers.txt:22` and is now false for the D ≥ 20 rows: those errors
belong to the old rows. The check still passes (72.63, 73.70) but its text says "~ 73".

## The remainder, recomputed (source: `logs_check_growth.txt:49-65`)

s = 2(D − 3). W²a² − [s + Ω₀ + Ω₁/s + …] times s^(k+1), from the csv rows. Ω₀ … Ω₃ are the closed forms of #26 and
#208 (R0: they equal `order5_output.txt`'s B4 values to 1e-50); Ω₄ = −12794.499860023 is #195's fit. ± is the rounding
of the 12th decimal given the crate's value, an enclosure; it is not the crate's own error (1e-60) and not a fit
spread.

| D | 3 terms × s² (old csv) | 3 terms × s² (new) | ± | 4 terms × s³ | ± | 5 terms × s⁴ | ± |
|---|---|---|---|---|---|---|---|
| 20 | 47.5456 | 47.5456 | 5.8e-10 | −912.0081 | 2.0e-08 | −8811.71 | 6.7e-07 |
| 50 | 66.0740 | 66.0740 | 4.4e-09 | −779.7619 | 4.2e-07 | −11930.65 | 3.9e-05 |
| 100 | 70.6708 | 70.6725 | 1.9e-08 | −717.1889 | 3.7e-06 | −12483.67 | 7.1e-04 |
| 200 | 72.6054 | 72.6308 | 7.8e-08 | −684.9892 | 3.1e-05 | −12666.76 | 1.2e-02 |
| 500 | 72.9870 | 73.6997 | 4.9e-07 | −665.6670 | 4.9e-04 | −12749.93 | 4.9e-01 |
| 1000 | 66.3152 | 74.0388 | 2.0e-06 | −659.2485 | 4.0e-03 | −12778.31 | 7.9e+00 |
| 2000 | −20.6375 | 74.2051 | 8.0e-06 | −656.0699 | 3.2e-02 | −12899.96 | 1.3e+02 |
| limit | | Ω₂ = 74.3694 | | Ω₃ = −652.8401 | | Ω₄ = −12794.4999 | |

- The three-term column tends to Ω₂ + Ω₃/s + Ω₄/s² + …, so it rises with a 1/D drift: 74.3694 − 0.3274 − 0.0032 =
  74.0388 at D = 1000.
- The four-term column tends to Ω₃ from above in modulus (−656.07 at D = 2000, −652.84 the limit; Ω₄/s = −3.2).
- The five-term column at D = 1000 and 2000 is limited by the 12-decimal rounding (± 7.9 and ± 130 against Ω₄ =
  −12794.5). Quote it for D ≤ 500, or use the 60-digit values of `order5_output.txt` (rem·s⁴ = −12773.08 at D = 1000,
  −12784.04 at D = 2000).
- R1 and R2 confirm the rows reproduce `w2_output.txt`'s and `order5_output.txt`'s printed 60-digit values within the
  enclosure (at D = 2000, R1: 0.029 against a bar of 0.032; R2: 1.2e+2 against 1.3e+2). At 1 decimal the D = 2000
  four-term value is ambiguous (−656.0 from the crate, −656.1 from the row); quote "−656".

## Every quoted number that reads growth.csv

### notes/wormhole_largeD/README.md

| line | old | new | source |
|---|---|---|---|
| 30 | "about 73/(2(D − 3))² at D = 200 and 500" | 72.63 (D = 200), 73.70 (D = 500) | `logs_paper_numbers.txt:18-19` |
| 31 | "D = 1000 gives 66 ± 6" | 74.04 | `logs_paper_numbers.txt:20` |
| 37-39 | rows "low by 2.7e-9 (D = 30) to 6e-6 (D = 2000)"; "66 ± 6 came from those low values" | rows replaced; the old "66 ± 6" is now 74.04 | `w2_output.txt:47-54`, above |
| 110-115 | W²a² here: 31.5751101276, 91.3457759034, 191.2830176842, 391.2543342170, 991.2379714266, 1991.2326553188 | 31.5751101284, 91.3457759136, 191.2830177304, 391.2543343808, 991.2379721479, 1991.2326572613 | `logs_paper_numbers.txt:14,16-20` |
| 119-120 | "values move by 4×10⁻¹⁰ to 1.5×10⁻⁶ ... (logs_modes_r.txt)" | true of `modes_r.py`, understates its error (low by 8e-10 to 1.9e-6, 6e-6 at D = 2000) | `w2_output.txt:47-54` |
| 121 | "fit in 1/D gives −4.77259 ± 10⁻⁵" | see below | raw `dev_06_fit_1overD.log` |

Lines 108, 109 (D = 4: 1.39708, D = 10: 12.0455) are the l-integration (`master.py modes`), not growth.csv: unchanged.
The TS13 and three-term columns are unchanged.

Proposed text, line 30-31:

```diff
-- Against the full numerics the three terms leave about 73/(2(D − 3))² at D = 200 and 500, as an O(1/D²) remainder
-  should. D = 1000 gives 66 ± 6, within that value's numerical error (table below).
+- Against #195's values (growth.csv, regenerated in #200) the remainder of the three terms times s² = (2(D − 3))² is 72.63,
+  73.70, 74.04 and 74.21 at D = 200, 500, 1000 and 2000. It tends to Ω₂ + Ω₃/s + … (74.04 = 74.37 − 0.33 − 0.003 at
+  D = 1000), as an O(1/D²) remainder should. The four-term remainder times s³ is −685.0, −665.7, −659.2 and −656 at the same D,
+  tending to Ω₃ = −652.84 (`logs_check_growth.txt`).
```

Proposed text, lines 37-39 (after "`w2_output.txt` has the 60-digit values."):

```diff
-  The `modes_r.py` values of the table below and of growth.csv (D ≥ 20) are low by 2.7e-9 (D = 30) to 6e-6
-  (D = 2000); `w2_output.txt` has the 60-digit values.
-  *Added by the coordinator:* the "66 ± 6" at D = 1000 above came from those low values.
+  The `modes_r.py` values were low by 8e-10 (D = 20) to 6e-6 (D = 2000). **Update from #200:** growth.csv's D ≥ 20 rows and the
+  table below now hold #195's quad-double values rounded to 12 decimals (`profiles.py`; `check_growth.py` checks them, `logs_check_growth.txt`).
+  The old "66 ± 6" at D = 1000 is 74.04.
```

Proposed table (lines 106-115) and the lines after it:

```diff
-| D | W²a² here | [TS13] | 2(D − 3) − 4 ln 2 + Ω₁/(2(D − 3)) |
+| D | W²a² (D ≤ 10: l-integration; D ≥ 20: #195's crate) | [TS13] | 2(D − 3) − 4 ln 2 + Ω₁/(2(D − 3)) |
-| 20 | 31.5751101276 | 31.5751101285 | 31.53398 |
-| 50 | 91.3457759034 | 91.3457759137 | 91.33830 |
-| 100 | 191.2830176842 | 191.2830177297 | 191.28114 |
-| 200 | 391.2543342170 | | 391.25387 |
-| 500 | 991.2379714266 | | 991.23790 |
-| 1000 | 1991.2326553188 | | 1991.23264 |
+| 20 | 31.5751101284 | 31.5751101285 | 31.53398 |
+| 50 | 91.3457759136 | 91.3457759137 | 91.33830 |
+| 100 | 191.2830177304 | 191.2830177297 | 191.28114 |
+| 200 | 391.2543343808 | | 391.25387 |
+| 500 | 991.2379721479 | | 991.23790 |
+| 1000 | 1991.2326572613 | | 1991.23264 |
-- The r-integration's values move by 4×10⁻¹⁰ to 1.5×10⁻⁶ from D = 20 to 1000 when its numerical parameters are changed
-  (`logs_modes_r.txt`).
+- `modes_r.py`'s values move by 4×10⁻¹⁰ to 1.5×10⁻⁶ from D = 20 to 1000 when its numerical parameters are changed
+  (`logs_modes_r.txt`). That understated its error: against #195's values it is low by 8e-10 to 1.9e-6 (6e-6 at D = 2000), and a
+  rerun on another machine moves its values by up to 3e-8 (#200).
```

The fit (line 121): no script in the repository produces "−4.77259 ± 10⁻⁵" (it came in with the first commit,
e9c2c47), and its design is not recorded. Polynomial fits of W²a² − 2(D − 2) in 1/D over D ≥ 50 (raw `dev_06`, not a
result): old rows −4.772593 (cubic), new rows −4.772588 (cubic) and −4.772589 (quartic), against −2 − 4 ln 2 =
−4.7725887. The printed "−4.77259 ± 10⁻⁵" stays true. Either keep it and say the fit is not scripted, or replace it by
a scripted number: the four-term series' relative difference from W²a² is 2.58e-12 at D = 2000
(`logs_check_growth.txt:75`).

README "Files" (line ~190-197): add `check_growth.py`, `logs_check_growth.txt`, and for `profiles.py` "growth.csv's
D ≥ 20 rows are #195's crate values".

### paper/3-wormhole-largeD/main.tex

| line | old | new | source |
|---|---|---|---|
| 138-143, col 2 | 31.5751101276, 91.3457759034, 191.2830176842, 391.2543342170, 991.2379714266, 1991.2326553188 | 31.5751101284, 91.3457759136, 191.2830177304, 391.2543343808, 991.2379721479, 1991.2326572613 | `logs_paper_numbers.txt:14,16-20` |
| 143, col 3 | 8.4×10⁻⁹ | 9.4×10⁻⁹ (9.352e-09) | `logs_check_growth.txt:64` |
| 142, col 4 | 73.0 | 73.7 | `logs_paper_numbers.txt:19` |
| 143, col 4 | 66 ± 6 | 74.0 | `logs_paper_numbers.txt:20` |
| 214 | "settles near 73 at D = 200-500", "66 ± 6" | see below | same |
| 116-117 | "moves W²a² by 4×10⁻¹⁰ at D = 20 and 1.5×10⁻⁶ at D = 1000" | true of `modes_r.py`; understates its error | `w2_output.txt:47-54` |

Other cells of the table do not change at the printed digits: column 3 at D = 20, 50, 100, 200, 500 (1.3×10⁻³,
8.2×10⁻⁵, 9.8×10⁻⁶, 1.2×10⁻⁶, 7.5×10⁻⁸) and column 4 at D = 20, 50, 100, 200 (48, 66, 71, 72.6). The parenthesized
TS13 values stay. The table now agrees with them to 4.8e-12, 9.4e-13 and 3.4e-12 (relative) at D = 20, 50, 100,
against 3.0e-11, 1.1e-10, 2.4e-10 before (`logs_check_growth.txt:38-40`). Abstract line 35 ("nine or more significant
figures", "up to D = 1000"): still true; the data run to D = 2000.

```diff
-20 & 31.5751101276 (31.5751101285) & $1.3\times10^{-3}$ & $48$ \\
-50 & 91.3457759034 (91.3457759137) & $8.2\times10^{-5}$ & $66$ \\
-100 & 191.2830176842 (191.2830177297) & $9.8\times10^{-6}$ & $71$ \\
-200 & 391.2543342170 & $1.2\times10^{-6}$ & $72.6$ \\
-500 & 991.2379714266 & $7.5\times10^{-8}$ & $73.0$ \\
-1000 & 1991.2326553188 & $8.4\times10^{-9}$ & $66\pm6$ \\
+20 & 31.5751101284 (31.5751101285) & $1.3\times10^{-3}$ & $48$ \\
+50 & 91.3457759136 (91.3457759137) & $8.2\times10^{-5}$ & $66$ \\
+100 & 191.2830177304 (191.2830177297) & $9.8\times10^{-6}$ & $71$ \\
+200 & 391.2543343808 & $1.2\times10^{-6}$ & $72.6$ \\
+500 & 991.2379721479 & $7.5\times10^{-8}$ & $73.7$ \\
+1000 & 1991.2326572613 & $9.4\times10^{-9}$ & $74.0$ \\
```

Caption (lines 123-126): "for $D\ge20$ they are from the $r$-integration of Sec. 2" no longer describes the rows. Proposed:

```diff
-$D\le11$ the values are Torii and Shinkai's~\cite{TS2013}; for $D\ge20$
-they are from the $r$-integration of Sec.~\ref{sec:numerics}, with Ref.~\cite{TS2013} in parentheses where available.
+$D\le11$ the values are Torii and Shinkai's~\cite{TS2013}; for $D\ge20$
+they are from a Taylor-series integration in quad-double arithmetic (the repository's \texttt{w2\_crate}), with Ref.~\cite{TS2013} in parentheses where available.
```

Text, lines 213-216:

```diff
-compare it with the numerical eigenvalues. The difference times $s^2$
-settles near $73$ at $D=200$--$500$, as an $O(D^{-2})$ remainder should; at $D=1000$ it is $66\pm6$, within that
-value's numerical uncertainty. The three terms are within $1.1\%$ ...
+compare it with the numerical eigenvalues. The difference times $s^2$
+is $72.6$, $73.7$, $74.0$ and $74.2$ at $D=200$, $500$, $1000$ and $2000$: it settles near $74$, with a slow $1/D$ drift, as
+an $O(D^{-2})$ remainder should. The three terms are within $1.1\%$ ...
```

If #210 adds the derived Ω₂, Ω₃: the remainder after the four-term series times s³ is −685.0, −665.7, −659.2, −656 at
D = 200, 500, 1000, 2000, tending to Ω₃ = −652.84 (table above); the five-term remainder times s⁴ is −12667 (D = 200),
−12750 (D = 500).

Sec. 2 (lines 116-120) and the method bullets (319-323): the sentence "Changing X from 40 to 50 ... moves W²a² by
4×10⁻¹⁰ at D = 20 and 1.5×10⁻⁶ at D = 1000" stays true of that integration, but its values are low by 8e-10 and 1.9e-6
against the crate, so it understates the error. "The second reproduces Torii and Shinkai's values to nine or more
significant figures" stays true of `modes_r.py` (worst 2.4e-10 relative at D = 100) and of the crate (worst 4.4e-10 at
D = 4, 4.8e-12 at D = 20). If the table's D ≥ 20 values are the crate's, "Two independent integrations" (line 319)
misdescribes the table's source: three computations exist. The wording is #210's.

Figure 1(a): `make_figures.py` `fig6_wormhole` reads the new csv without change (run into a scratch folder, raw
`run_05`; `source.startswith("TS13")` still splits the markers). The largest ordinate change is +5.95e-6 (D = 2000) on
an axis of 2.6, so nothing visible moves. The PDF was not regenerated. The caption "this work's for D = 20-2000" can
say "#195's".

### Other places

| file:line | what | old | new / action |
|---|---|---|---|
| `one_mode.py:145` | hard-coded W² for N4 | 31.5751101276, 91.3457759034 | 31.5751101284, 91.3457759136. N4's bar is 1e-6, so its verdict cannot change (old and new differ by 2.5e-11 and 1.1e-10 relative). Not run. |
| `one_mode.log:33,35` | "against W^2 = 31.5751101276 (relative 2.4e-11)", "= 91.3457759034 (relative 1.1e-10)" | | would be rerun after the edit above |
| `one_mode.md:125` | "to 10⁻⁶" | | unchanged |
| `notes/wormhole_next.md:241` | "The remainder is about 73/(2(D − 3))²" | | dated planning text; superseded by README's #26 update (74.37 and the drift above) |
| `omega2.md:169` | "growth.csv is not used (#195: its large-D rows are low by 8e-10 to 6e-6)" | | true as history; append "(regenerated in #200)" |
| `w2_highprec.md:120-139, 221, 242-245, 300`, `w2_output.txt:37-54,171-176` | the old rows beside the crate (G1, tight run) | | records of the run at 514ac42; leave. A driver rerun now reads the new csv: G1 gives 8 of 16 (the 8 TS13 rows still fail, #202) and the "growth.csv (modes_r.py)" column and "tight − csv" become meaningless. |
| `order4.py:808` | `read_w2` regex needs the label `modes_r.py` or `TS13 Table I` in `w2_output.txt`'s table | | breaks on a `w2_output.txt` regenerated after this change (the label becomes `w2_crate qd (#195)`). `profiles.py` and `check_growth.py` read the table under any label. |

## Reported, not fixed

- **The comma.** The task's source string broke `w2_highprec.py`'s reader; the coordinator ruled the comma-free string
  (above).
- **`profiles.csv`.** Not rewritten in the commit. `profiles.py` still takes the profiles' eigenvalue at D = 20, 100,
  1000 from `modes_r.py`, which is not growth.csv's. Using the crate's would change the profiles by at most 9.5e-10,
  1.7e-9 and 7.7e-9 (relative; raw `dev_05`), below the figure's resolution; the rerun on rpi5 alone differs from the
  committed file by up to 1.8e-9. If the profiles should carry the table's value, replace `m.W2_of(D)` by
  `float(round12(crate_w2()[D]))` in `write_profiles`. Decide.
- **TS13 rows (D = 4 … 11).** Unchanged. They differ from the crate by +6.19e-10 (D = 4) to −2.46e-10 (D = 10); #202
  decides their digits.
- **The 1/D fit in README.md line 121** has no script (above).
- **R2 at D = 2000** passes with 1.2e+2 against a bar of 1.3e+2: the csv's rounding alone makes the five-term column
  there a one-digit number.

## Checks and controls

Written in the raw README before the first run; targets declared from a development run (`dev_04`) before the runs of
record (`run_*`).

| check | what | result |
|---|---|---|
| G1 | header and the 16 rows | pass |
| G2 | 8 rows D ≥ 20 equal the crate rounded to 12 decimals | 8 of 8 |
| G3 | source strings | pass |
| G4 | 8 TS13 rows equal profiles.py's TS13 values at 12 decimals | 8 of 8 |
| G5 | table and list agree at D = 50 | pass |
| G6 | TS13 rows equal 613e866's | 8 of 8 |
| R0 | closed forms Ω₀ … Ω₃ vs `order5_output.txt` | differences 8.3e-109 to 4.6e-59 (bar 1e-50) |
| R1, R2 | remainder columns vs the printed 60-digit values | pass (above) |

56 mutations (`last_digit`, `drop_row`, `source` at all 16 D; `ts13_swap`, `crate_scale`, `round_trunc`, `list50`,
`wrong_om2`, `wrong_om3`, `s_off`): all on their declared targets in `run_03`, `run_07` and `run_09` (the last for the
final code), none silent. Shared by the checker and `profiles.py`: the input `w2_output.txt` and TS13's values. The
parsers and the roundings (exact rationals against `Decimal`) are separate. Rerun: `python3 profiles.py
--growth-only`, `python3 check_growth.py > logs_check_growth.txt`, `python3 paper_numbers.py >
logs_paper_numbers.txt`, from `notes/wormhole_largeD`.

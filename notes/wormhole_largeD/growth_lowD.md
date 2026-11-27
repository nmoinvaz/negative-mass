# #221: growth.csv's D = 4 … 11 rows from #195's crate

Status: **checked** by `check_growth.py` (G1 to G6, R0 to R2 pass; 64 mutations and the unmutated run, 65 of 65 on their
declared targets). All 16 rows of growth.csv are now #195's quad-double values rounded to 12 decimals, source
`w2_crate qd (#195)`. **None of the 31 printed numbers in Paper III and README.md that read the eight rows changes.**
`ts13_digits.py` was added to the file list afterwards (hc, 2026-11-27): its C1 failed on the new csv, so C1 now reads
the old rows as a constant and a new C2 checks the current rows against the crate; its `--selftest` passes (below).
`make_figures.py`'s marker split and the text proposals below go into #210. Branch `wip/i221-growth-lowD`, raw
`negative-mass-raw/I221-rpi5` (README first).

## Result

| acceptance | verdict |
|---|---|
| `check_growth.py` checks all 16 rows against the crate, and its row mutation fires on each | **met**: G2 passes 16 of 16; `last_digit:D` fires G2 at each of the 16 D (`logs_check_growth.txt:3`, raw `run_03_mutations/`) |
| `paper_numbers.py` reruns with "checks failed: 0", listing every changed number | **met**: exit 0, `checks failed: 0` (`logs_paper_numbers.txt:37`). 8 of 37 lines changed: the W²a² column at 10 decimals for D = 4, 5, 6, 7, 9, 10, 11 and the source column for all 8 (below) |
| The runner expects none of Paper III's printed numbers to change; any that does is listed | **met**: 0 of 31 printed cells change (below). The caption and the Fig. 1 caption name TS13 as the source of the D ≤ 11 values, which is no longer true: text, not numbers |
| growth.csv keeps three columns | **met**: G1; `w2_highprec.py`'s `read_csv` reads the file (16 rows, raw `run_08_read_csv.log`) |

## What was done

```
notes/wormhole_largeD/
├── profiles.py            write_growth(): all 16 rows from crate_w2(), round12() half-even, source "w2_crate qd (#195)";
│                          the TS13 dict stays (write_profiles reads it for the D = 10 profile; ts13_digits.py T2 and
│                          check_growth.py G4 read its text); profiles.csv not rewritten
├── growth.csv             lines 2 to 9 changed (D = 4 … 11); lines 10 to 17 (D >= 20) byte-identical to 4bced01
├── check_growth.py        G2, G3 over all 16 rows; G4 = D <= 11 rows within relative 1e-8 of TS13's printed values;
│                          G6 = D >= 20 rows equal 4bced01's; enclosure 5e-13 s^p for all rows; mutations revert_lowD,
│                          ts13_digit:D (new), source:D relabelled
├── logs_check_growth.txt  its output
├── logs_paper_numbers.txt paper_numbers.py rerun
├── ts13_digits.py         C1 against OLD_ROWS (the old csv rows, a constant), new C2 (current rows against the crate),
│                          mutations csv retargeted, old_rows and crate_lowD new
├── ts13_digits_output.txt its output and `--selftest`, refreshed
└── growth_lowD.md         this note
```

## The eight rows, old and new (source: `growth.csv:2-9`, `git show 4bced01:notes/wormhole_largeD/growth.csv`)

| D | old csv (TS13 rounded) | new csv (crate) | new − old | relative |
|---|---|---|---|---|
| 4 | 1.397052433715 | 1.397052434334 | +6.19e-10 | 4.4e-10 |
| 5 | 2.984958930278 | 2.984958930352 | +7.40e-11 | 2.5e-11 |
| 6 | 4.686620542995 | 4.686620543119 | +1.24e-10 | 2.6e-11 |
| 7 | 6.462584141263 | 6.462584141370 | +1.07e-10 | 1.7e-11 |
| 8 | 8.289759363063 | 8.289759363072 | +9.00e-12 | 1.1e-12 |
| 9 | 10.153553045187 | 10.153553044946 | −2.41e-10 | 2.4e-11 |
| 10 | 12.044265014744 | 12.044265014498 | −2.46e-10 | 2.0e-11 |
| 11 | 13.955209167665 | 13.955209167496 | −1.69e-10 | 1.2e-11 |

Each new row is within 4.93e-13 of the crate's value (`logs_check_growth.txt:11-18`; half a unit of the 12th decimal is 5e-13);
the closest a crate value comes to a rounding tie is 0.007 of a unit (D = 4), far above the crate's own error. Against TS13's
printed digits the new rows differ by 4.4e-10 (D = 4) or less, relative (G4, `logs_check_growth.txt:5`). The old rows were
TS13's printed values to within 4.1e-13 (`logs_check_growth.txt:30-37`).

## Numbers that read the eight rows

`printed_digits.py` (raw folder) parses each printed number from `main.tex` and README.md, recomputes it from the old and the new
rows at the printed digits, and compares. 31 cells, 31 unchanged (raw `run_09_printed_digits.log`); its mutation `bump:D` fires
exactly the declared cells, 9 of 9 runs on target (raw `run_10_compare.log`).

| where | cells | change at the printed digits |
|---|---|---|
| `paper/3-wormhole-largeD/main.tex:132-137`, Table 1, D = 4, 5, 6, 7, 8, 10 | W²a² (8 decimals), relative difference from the three-term formula, (difference)·s²: 18 | none |
| `main.tex:280-283`, Table 2, D = 4, 5, 6 | 1/(Wa) and the predicted slope, 3 decimals: 6 | none |
| `notes/wormhole_largeD/README.md:68-69, 76-77`, e-folding times 0.846, 0.579, 0.462 | 5 (the 0.462 of lines 76-77 is split across the lines; it is line 69's) | none |
| `README.md:110-111`, [TS13] column at D = 4 (1.39705243) and D = 10 (12.0442650) | 2 | none |

The same holds for the rest of `paper_numbers.py`'s output: D ≥ 20, the e-folding table, the ST17 slopes and sigmas, the SH2002
comparison, and its own checks (none fails). D = 9 is printed nowhere: at 8 decimals it would change, 10.15355305 to
10.15355304. D = 11 is printed nowhere and does not change at 8 decimals (13.95520917). README.md's column "W²a² here" at D = 4
and 10 (1.39708, 12.0455) is the l-integration, not growth.csv. The abstract's "nine or more significant figures" (`main.tex:35`,
119-120) stays true: the crate's rows agree with TS13's printed digits to 9.35 to 12.03 digits (`ts13_digits.py` F1).

`paper_numbers.py` lines that changed (`logs_paper_numbers.txt:6-13`; the remaining columns are equal at the printed digits):

| D | W²a² old | W²a² new | source column |
|---|---|---|---|
| 4 | 1.3970524337 | 1.3970524343 | `TS13 Table I` → `w2_crate qd (#195)` |
| 5 | 2.9849589303 | 2.9849589304 | same |
| 6 | 4.6866205430 | 4.6866205431 | same |
| 7 | 6.4625841413 | 6.4625841414 | same |
| 8 | 8.2897593631 | 8.2897593631 | same (the value is equal at 10 decimals) |
| 9 | 10.1535530452 | 10.1535530449 | same |
| 10 | 12.0442650147 | 12.0442650145 | same |
| 11 | 13.9552091677 | 13.9552091675 | same |

## Text that names the source and is now wrong (proposed)

Paper III, `main.tex:124-126`, Table 1's caption. This also fixes the D ≥ 20 half that #200's note proposed and #210 left:

```diff
-\caption{$W^2a^2$ for the growing mode. For $D\le11$ the values are Torii and Shinkai's~\cite{TS2013}; for $D\ge20$
-they are from the $r$-integration of Sec.~\ref{sec:numerics}, with Ref.~\cite{TS2013} in parentheses where available.
+\caption{$W^2a^2$ for the growing mode, from a Taylor-series integration in quad-double arithmetic (the repository's
+\texttt{w2\_crate}). Ref.~\cite{TS2013}'s printed values are in parentheses for $D\ge20$; for $D\le11$ they equal these to the
+eight decimals shown.
```

"Equal to the eight decimals shown": TS13's printed values rounded to 8 decimals equal the new rows at D = 4, 5, 6, 7, 8, 10
(and 11; `ts13_digits_output.txt:54-62`, column "A = B at 8 decimals"); D = 9 would not.

Paper III, `main.tex:225-227`, Fig. 1 caption, if `make_figures.py` keeps the open circles for D ≤ 11 (below):

```diff
-(a) $W^2a^2-2(D-3)$ against $1/(D-3)$: Torii and Shinkai's eigenvalues~\cite{TS2013} for $D=4$--$11$ (open
-circles) and this work's for $D=20$--$2000$ (filled), with the two-term (dashed) and three-term (solid) large-$D$
+(a) $W^2a^2-2(D-3)$ against $1/(D-3)$: the eigenvalues for $D=4$--$11$ (open circles; they agree with Torii and
+Shinkai's~\cite{TS2013} to $4.4\times10^{-10}$) and $D=20$--$2000$ (filled), with the two-term (dashed) and three-term (solid) large-$D$
```

`paper/figures/make_figures.py:247` splits the markers on `source.startswith("TS13")`. No source starts with it now, so every
marker takes the "this work" style and the legend keeps an entry for an empty series. Run into a scratch folder: it runs; the
committed and the new PNG are in raw (`run_07_fig6_committed.png`, `run_07_fig6_new.png`). The ordinates move by at most 6.2e-10
on an axis of 2.6, so no point moves. Not edited (not on the file list); the PDF in the paper is not regenerated. Proposal,
which restores the committed figure and keeps the two marker styles by dimension:

```diff
-    ts = np.array([r["source"].startswith("TS13") for r in rows])
+    ts = D <= 11
-    ax.plot(x[ts], y[ts], "o", ms=4.5, mfc="white", mec=ORANGE, mew=1.2, label="Torii & Shinkai (2013)")
-    ax.plot(x[~ts], y[~ts], "o", ms=3.5, color=BLUE, label="this work")
+    ax.plot(x[ts], y[ts], "o", ms=4.5, mfc="white", mec=ORANGE, mew=1.2, label=r"$D\le11$")
+    ax.plot(x[~ts], y[~ts], "o", ms=3.5, color=BLUE, label=r"$D\ge20$")
```

`notes/wormhole_largeD/README.md`, "Files" (line 202), after #200's proposal:

```diff
 - `profiles.py`: the mode profiles and the table of growth rates behind Paper III's figure (`profiles.csv`,
-  `growth.csv`).
+  `growth.csv`: all 16 rows are #195's quad-double values rounded to 12 decimals, #221).
+- `check_growth.py`: checks `growth.csv` against `w2_output.txt`, and the remainder columns; `--mutate=NAME`
+  (`logs_check_growth.txt`, #200, #221).
```

`notes/literature.md:94-96`, the TS13 entry:

```diff
-    od lanes agree to 60 digits but share code). growth.csv's D = 4 … 11 rows are Table I rounded to 12 decimals, so they
-    carry that error; no csv digit there goes beyond the printed ones. README.md's quoted [TS13] values are roundings of
-    Table I.
+    od lanes agree to 60 digits but share code). growth.csv's D = 4 … 11 rows were Table I rounded to 12 decimals (#202); since
+    #221 they are the crate's values, and Table I's printed values are `ts13_digits.py`'s cross-check. README.md's quoted
+    [TS13] values are roundings of Table I.
```

`notes/wormhole_largeD/w2_highprec.md`, after the "Update from #202" paragraph (line 141 ff.), the csv is no longer what the
paragraph describes:

```diff
+*Update from #221:* growth.csv's D = 4 … 11 rows now hold this crate's values (`w2_crate qd (#195)`). A driver rerun reads the
+new csv: G1 is 16 of 16 within 5e-13 (the largest |crate − csv| is 4.93e-13 at D = 4, `logs_check_growth.txt:11-26`; the driver was not
+rerun), and the table's source column reads `w2_crate qd (#195)` for all 16 rows.
```

`notes/wormhole_largeD/paper_numbers.py`, docstring (supersedes #200's proposal for line 5; its other two stale lines stay:
the print at `logs_paper_numbers.txt:22` and the check text "~ 73"):

```diff
-Inputs: growth.csv (profiles.py: Torii & Shinkai's Table I for D <= 11, modes_r.py for D >= 20); logs_order3.txt
+Inputs: growth.csv (profiles.py: #195's quad-double values from w2_output.txt, rounded to 12 decimals); logs_order3.txt
```

## Reported, not fixed

- **`ts13_digits.py` (added to the file list by the coordinator after the first report).** On the new csv its C1 ("the csv's D = 4 … 11
  rows equal −ω² at 12 decimals") failed 0 of 8, and `--selftest` failed at its baseline (raw `run_05_*`). Now: C1 checks TS13's
  printed values, rounded to 12 decimals, against `OLD_ROWS`, a constant holding growth.csv's rows at 4bced01 (equal to the git file at
  4bced01 and 613e866, raw `dev_06_verify_old_rows.log`); it does not read the current csv. C2 checks the current rows D = 4 … 11
  against the crate rounded to 12 decimals, as `check_growth.py` G2 does, with the script's own parser and `Decimal` rounding. The
  baseline passes (C1 8 of 8, C2 8 of 8, the other nine checks as before) and `--selftest` passes, 18 of 18 runs (raw `run_12_*`,
  `ts13_digits_output.txt`). One declared target changed, `csv` (C1 to C2); `old_rows` (C1) and `crate_lowD` (C2) are new; the other
  fifteen are as in #202. Predictions, development run and runs of record: raw README, amendments A3 and A3b. `README.md`'s
  table is still read by R1 (the [TS13] column, `read_readme`), so #210's edit of that table must keep its layout.
- **`make_figures.py`'s marker split**, above.
- **Readers of `w2_output.txt` keyed on the source label.** `w2_output.txt` is unchanged, so nothing breaks now. A rerun of
  `w2_highprec.py` prints `w2_crate qd (#195)` in the table's source column. Then `notes/black_universe/stability.py:54`
  (`planted_value`, regex on the row labelled `TS13 Table I` at D = 4) raises `RuntimeError`, and `order4.py:808` matches 0 of 16
  rows (it already missed the D ≥ 20 rows after #200). Tested on a relabelled copy (raw `run_11_relabel_readers.log`, scratch copy not kept). `check_growth.py`
  and `profiles.py` read the table under any label.
- **Prediction miss in the plan.** I predicted `source:D` at D ≥ 20 fails G3 alone; it fails G3 and G6, since G6 compares the whole
  row (raw README, amendment A1). No criterion or bar was changed after a run.
- **Two records of TS13's digits.** `profiles.py`'s `TS13` dict (a transcription; `write_profiles` uses it for the D = 10 profile) and
  `ts13_digits.py`'s `PRINTED` (typed from the TeX source) are equal (T2). The dict is kept because three readers use it; growth.csv does not.
- **`profiles.csv`** is unchanged and `profiles.py` was run with `--growth-only`. Its D = 10 profile still uses TS13's eigenvalue
  (12.0442650147438; the crate's is 12.044265014498, 2.5e-10 lower) and D = 20, 100, 1000 use `modes_r.py`'s (#200's open decision).
  Neither shows at the figure's resolution.
- **Other prose that reads TS13's or the crate's D ≤ 11 values at 9 digits or fewer** (`notes/critical_ghost.md:93,103`,
  `notes/ellis_efold.md:21`, `notes/black_universe/README.md:60,99-100`) quotes TS13 or `w2_output.txt`, not growth.csv: unchanged.

## Checks and controls

Written in the raw README before the first run of the new checker (criteria), predictions before the development run, targets
declared from it (`dev_03`) before the runs of record (`run_03`). Shared by `check_growth.py` and `profiles.py`: the input
`w2_output.txt` and TS13's dict; the parsers and the roundings (exact rationals against `Decimal`) are separate. The csv rows are
also supported by code that shares nothing with the crate: `one_mode.log` equals the crate to 12 digits at D = 4, 10, 20, 50
(`ts13_digits.py` O1, passes).

| check | what | result |
|---|---|---|
| G1 | header and the 16 rows, three fields each | pass |
| G2 | all 16 rows equal the crate's qd value rounded to 12 decimals | 16 of 16 |
| G3 | source `w2_crate qd (#195)` on all 16 rows | pass |
| G4 | D ≤ 11 rows within relative 1e-8 of TS13's printed values (K1's bar) | 8 of 8, worst 4.43e-10 |
| G5 | table and list agree at D = 50 | pass |
| G6 | D ≥ 20 rows equal 4bced01's | 8 of 8 |
| R0 | closed forms Ω₀ … Ω₃ against `order5_output.txt` | 8.3e-109 to 4.6e-59 (bar 1e-50) |
| R1, R2 | remainder columns against the printed 60-digit values | pass |
| `ts13_digits.py` C1 | TS13's printed values rounded to 12 decimals equal `OLD_ROWS` (the old csv rows), which print 13 to 14 decimals | 8 of 8 |
| `ts13_digits.py` C2 | the current csv's D = 4 … 11 rows equal the crate rounded to 12 decimals | 8 of 8 |

| controls | result |
|---|---|
| `check_growth.py --mutate=`: `last_digit:D`, `drop_row:D`, `source:D` at all 16 D; `ts13_digit:D` at D = 4 … 11; `revert_lowD`, `ts13_swap`, `crate_scale`, `round_trunc`, `list50`, `wrong_om2`, `wrong_om3`, `s_off` | 64 of 64 fire exactly their declared checks, none silent (`run_03_compare_targets.log`) |
| `ts13_digits.py --mutate=`: 15 of #202 (`csv` now targets C2) and `old_rows`, `crate_lowD` (new) | 17 of 17 fire exactly their declared checks, and `--selftest` passes (raw `dev_05_*`, `run_13_*`) |
| `printed_digits.py --mutate=bump:D`, D = 4 … 11 | 6 fire the declared cells, `bump:9` and `bump:11` fire nothing by declaration (D printed nowhere); 9 of 9 on target |

`revert_lowD` is the change not applied (the eight rows back to TS13's rounded values and label): it fires G2 and G3, and G4 does
not, since TS13's values are within its bar. `ts13_digit:D` is what G4 alone catches. Rerun from `notes/wormhole_largeD`:
`python3 profiles.py --growth-only`, `python3 check_growth.py > logs_check_growth.txt`,
`python3 paper_numbers.py > logs_paper_numbers.txt`.

## What ran

| command | exit |
|---|---|
| `profiles.py --growth-only` (twice; the second left growth.csv byte-identical) | 0 |
| `check_growth.py` (unmutated) | 0 |
| `check_growth.py --mutate=…` × 64 | 1 each, on target |
| `paper_numbers.py` | 0 |
| `printed_digits.py` (raw) and its 8 mutations | 0; 1 for the six that fire, 0 for `bump:9`, `bump:11` |
| `ts13_digits.py` (retargeted), `--selftest` (18 runs), and its 17 `--mutate=` runs | 0, 0; 1 for each of the 17 mutations, 18 of 18 on target (raw `run_12_*`, `run_13_*`) |
| `make_figures.fig6_wormhole` into a scratch folder; `w2_highprec.read_csv` | 0, 0 |

## Recommended tasks

None open from this issue. By the coordinator's ruling, `make_figures.py` (the marker split above) and every paper or README text
(the diffs above) go into #210. `stability.py:54` and `order4.py:808` stay as they are; the dependency on `w2_output.txt`'s source
label is noted under "Reported, not fixed".

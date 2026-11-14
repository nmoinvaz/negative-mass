# References check (issue #49)

Run on 2026-11-13 on dh610 for [#49](https://github.com/nmoinvaz/negative-mass/issues/49). Every bibliographic entry of
Papers I–III and every journal detail of `notes/novelty/RESULTS.md`, `notes/critical_ghost.md`, `notes/bh_to_wormhole.md`
and `notes/wormhole_largeD/one_mode.md` was compared with arXiv, INSPIRE and the publisher's DOI record. The tables at
the end list every entry with the source used. The three papers and four notes hold 142 table rows, about 115
distinct works, since the papers and notes share many references.

## Method

- **arXiv.** The API record of each id (title, authors, journal-ref, DOI, version), and the abstract pages of
  2608.24577 v1 and v2.
- **INSPIRE.** The literature API by arXiv id (publication info, DOIs, titles, authors). The record number is in the
  tables.
- **Publisher.** The Crossref record of the DOI, which is the publisher's own deposit, for the journal, volume, pages,
  year, title and authors. The Oxford University Press page was read for Gibbons–Hartnoll–Ishibashi, whose title
  differs from arXiv's.
- **Books and items without a DOI.** Crossref, Open Library and zbMATH.
- **Comparison.** For every journal citation a script compared volume, first page or article number, and year with
  INSPIRE and with Crossref, and a citation passes when at least one of the two agrees. Titles and author lists were
  compared by script (normalised for case, hyphens, dashes and accents) and by eye. Not counted as mismatches: case,
  ". I." against ": I.", hyphens against en dashes in compound names, and ASCII transliterations of accents.
- **The check can fail.** The same script with every volume number raised by one flags all 104 citations that carry a
  volume. On the entries as corrected it flags none. It found no mismatch of volume, page or year in the original
  entries, so every correction below is a title, an author, a missing journal reference or an unverified mark. The
  title comparison flagged GGS09c and JJN23, and Kain2021 was found by comparing with INSPIRE's title.
- **Not read.** Only bibliographic details were checked. The content that a note quotes from a source was not re-read.
- **Raw data.** The API responses are not stored in the repository. They are reproducible from the ids and DOIs in the
  tables.

## Result

| File | Rows | Corrected | *(unverified)* marks resolved | Still unverified |
|---|---|---|---|---|
| paper/1-thresholds/main.tex | 17 | 6 | 0 | 0 |
| paper/2-boson-star/main.tex | 19 | 3 | 0 | 0 |
| paper/3-wormhole-largeD/main.tex | 12 | 0 | 0 | 0 |
| notes/novelty/RESULTS.md | 36 | 3 | 6 (the "Limits of this search" bullet) | 0 |
| notes/critical_ghost.md | 39 | 3 | 10 | 0 |
| notes/bh_to_wormhole.md | 13 | 2 | 1 | 0 |
| notes/wormhole_largeD/one_mode.md | 6 | 0 | 3 | 0 |

The two companion-paper entries (Moinvaziri2026a and Moinvaziri2026b) are this repository's own papers. They have no
external record, and their titles match the `\title` of the paper they name. No entry is left *(unverified)*.

## Corrections

Every mismatch found. A mark that was only confirmed is not listed here, it is in the tables.

| # | File | Was | Now | Source |
|---|---|---|---|---|
| 1 | Papers I, II, Shirokov2026 | `N. M. Shirokov`, title "positive--negative" | `N. M. Shirokov and I. Nachevsky`, "positive-negative" | arXiv abs v1 (one author) and v2 (two authors), INSPIRE 3195733 |
| 2 | Papers I, II, Alves2025 | `M. de Aguiar Alves` | `N. de Aguiar Alves` (Níckolas) | arXiv, INSPIRE 2813778, Crossref |
| 3 | Paper I, GHI2005 | "On the stability of naked singularities" | "On the stability of naked singularities with negative mass" | Crossref, OUP page (arXiv and INSPIRE have the shorter title) |
| 4 | Paper I, MP2014 | "de Sitter space-time" | "de Sitter spacetime" | Crossref (APS title) |
| 5 | Paper I, Johnson2020 | arXiv only | adds Gen. Relativ. Gravit. 52, 80 (2020) | INSPIRE 1757530, Crossref 10.1007/s10714-020-02732-9 |
| 6 | Paper I, Zerilli1970 | one title for PRL 24, 737 and PRD 2, 2141 | adds the PRD title "Gravitational field of a particle falling in a Schwarzschild geometry analyzed in tensor harmonics" | Crossref for both DOIs |
| 7 | Paper II, Kain2021 | "Radial oscillations and stability of multiple-field boson stars" | "Boson stars and their radial oscillations" | INSPIRE 1866961, Crossref (the old title is that of Kain's PRD 102, 023001 (2020) on multiple-fluid compact stars) |
| 8 | RESULTS.md, GHI | `PTP 113, 963` | adds `(2005)` | Crossref |
| 9 | RESULTS.md, Johnson | `GRG 52` | `GRG 52, 80` | INSPIRE, Crossref |
| 10 | RESULTS.md, Tamaki–Sakai | arXiv id only, label 2011 | adds `PRD 81, 124041 (2010)` | INSPIRE 860321, Crossref (arXiv posted 2011, paper published 2010) |
| 11 | RESULTS.md, limits bullet | six journal references "not checked" | states they were checked and correct | this file |
| 12 | critical_ghost.md, GGS09c | "On the instability of charged wormholes ..." | "Instability of charged wormholes ..." | Crossref 10.1103/PhysRevD.80.024023 |
| 13 | critical_ghost.md, JJN23 | `Jiménez-Vázquez`, "Confinement of exotic matter. I. Static solutions" | `E. Jiménez-Vázquez`, "Confinement of exotic matter: static solutions" | Crossref 10.1103/PhysRevD.107.064061 |
| 14 | critical_ghost.md, DMW23 | `P.-B. Ding, T.-X. Ma, Y.-Q. Wang` | `P.-B. Ding, T.-X. Ma, T.-F. Fang, Y.-Q. Wang` | Crossref 10.1007/JHEP04(2024)033, INSPIRE 2664190 (arXiv v1 has three authors) |
| 15 | critical_ghost.md, "Journal references" bullet | "marked *(unverified)* were not checked" | states they were checked | this file |
| 16 | bh_to_wormhole.md, HKL01 | `(2D dilaton gravity; journal *(unverified)*)` | `PRD 65, 064003 (2002)` | INSPIRE 564565, Crossref |
| 17 | bh_to_wormhole.md, Ter22 | `(journal not found)` | `PRD 106, 044035 (2022)` | INSPIRE 2048052, Crossref |

The marks resolved without a mismatch were GGS09c (PRD 80, 024023), DHNS09 (IJMPD 18, 1665), ST17 (PRD 96, 044009),
BFZ11 (EPJC 71, 1791), BKZ12 (PRD 86, 024028), DFKKR13 (PRD 88, 124028), Fro97 (PRD 56, 6433), HK04 (PRD 70, 101502),
Ber26 (PRD 113, 104041), DMW23 (JHEP 04 (2024) 033), the six references of the RESULTS.md bullet, and Wei87, BS91 and
Sim05 of `one_mode.md`. The sourcing remarks "(journal from citations)" and "(journal from [H09]'s reference list)" in
`bh_to_wormhole.md` are removed, since both journals are confirmed directly.

## What the sources could not settle

- **Moinvaziri2026a, Moinvaziri2026b.** Unpublished, so nothing external exists. Sought: arXiv and INSPIRE, no record.
- **Wei87.** Springer's own page does not show the series volume to a non-browser fetch. The volume 1258 comes from
  zbMATH, and the publisher, series and year from Crossref.
- **ReedSimon2, BirmanSolomyak, BS91.** Checked in library catalogues (Open Library, Crossref's Springer deposit), not
  on a publisher page. The place "Dordrecht" for BS91 is from the note and was not checked separately.
- **H02.** A conference contribution. INSPIRE lists it as a conference paper with no journal, so "journal not found"
  stands.

## Outside the files this task may edit

For the coordinator. Nothing below was edited.

1. **notes/warp/review.md, lines 138, 247, 248 (Aubin, Talenti).** Both verified, so the *(unverified)* marks can go.
   - Aubin, "Problèmes isopérimétriques et espaces de Sobolev", J. Differential Geom. 11(4), 573–598 (1976), Crossref
     10.4310/jdg/1214433725.
   - Talenti, "Best constant in Sobolev inequality", Ann. Mat. Pura Appl. 110, 353–372 (1976), Crossref
     10.1007/BF02418013.
   - Only the citation was checked. That the sharp constant is S = 3(π/2)^{4/3} and is attained by (1 + r²)^{−1/2} was
     not re-read in either paper.
2. **FL22 (notes/warp/mass_bounds.md:114, prior_work.md:199, literature.md:147, frank_loss.md:301).** Crossref
   10.4171/jems/1475: J. Eur. Math. Soc. 28, 305–331. Volume and pages match the notes. Crossref dates the article 2024
   (online first), so the year 2026 is not confirmed by that record.
3. **Other warp and wormhole notes.** Martín-Moruno–Visser 1802.00865 is CQG 35, 125003 (2018) (INSPIRE 1653135).
   Frauendiener–Nester–Szabados 1105.5008 is CQG 28, 185004 (2011). `notes/wormhole_next.md:564` leaves the journal of
   [NM08] *(unverified)*, and it is CQG 25, 055009 (2008).
4. **Shirokov–Nachevsky v2.** arXiv:2608.24577 v2 (2026-10-26) has a different abstract from v1 (2026-08-25). v1 reports
   "thirty-four evolutions", a speed of 0.056c by t = 400 and a lone star surviving to t = 1000. v2 reports a pair that
   reaches a tenth of the speed of light. The numbers that Paper II and `notes/novelty/RESULTS.md` take from Shirokov,
   such as the slow drift of a lone phantom star and the mass ladder, were read in v1. See Recommended tasks on the
   issue.
5. **PDFs.** `paper/1-thresholds/main.pdf` and `paper/2-boson-star/main.pdf` are not rebuilt. dh610 has no LaTeX, so the
   `.tex` edits were not compiled. They change only text inside `\bibitem` entries.
6. **Published since.** arXiv-only entries whose journal now exists: Kain26 (CQG 43, 045014 (2026)), Caridi–Corelli–Pani
   (PRD, DOI 10.1103/3c7y-437m), Trivedi–Loeb (Phys. Dark Univ. 53, 102409 (2026)), Raulot 2509.05244 (CQG 43, 045011
   (2026)). The notes do not claim a journal for them, so none is wrong.

## Entries checked

### paper/1-thresholds/main.tex (Paper I)

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| Bondi1957 | H. Bondi, Rev. Mod. Phys. 29, 423 (1957) | Crossref 10.1103/RevModPhys.29.423 | verified (volume, page and year match in Crossref) |
| Shirokov2026 | N. M. Shirokov, I. Nachevsky, arXiv:2608.24577 | arXiv:2608.24577v2; INSPIRE 3195733; arxiv.org/abs/2608.24577v1 and v2 | corrected, authors and title hyphen (v1 of 2026-08-25 is by Shirokov alone, v2 of 2026-10-26 adds Nachevsky and writes "positive-negative") |
| GHI2005 | Gibbons, Hartnoll, Ishibashi, Prog. Theor. Phys. 113, 963 (2005), hep-th/0409307 | arXiv:hep-th/0409307v1; INSPIRE 660779; Crossref 10.1143/PTP.113.963; academic.oup.com/ptp/article-lookup/doi/10.1143/PTP.113.963 | corrected, title, the published title is "On the stability of naked singularities with negative mass" (arXiv and INSPIRE carry the shorter one) |
| GD2006 | Gleiser, Dotti, Class. Quantum Grav. 23, 5063 (2006), gr-qc/0604021 | arXiv:gr-qc/0604021v2; INSPIRE 713860; Crossref 10.1088/0264-9381/23/15/021 | verified (volume, page and year match in INSPIRE and Crossref) |
| DG2009 | Dotti, Gleiser, Class. Quantum Grav. 26, 215002 (2009), 0809.3615 | arXiv:0809.3615v3; INSPIRE 797174; Crossref 10.1088/0264-9381/26/21/215002 | verified (volume, page and year match in INSPIRE and Crossref) |
| CC2006 | Cardoso, Cavaglià, Phys. Rev. D 74, 024027 (2006), gr-qc/0604101 | arXiv:gr-qc/0604101v2; INSPIRE 715170; Crossref 10.1103/PhysRevD.74.024027 | verified (volume, page and year match in INSPIRE and Crossref) |
| Dotti2022 | G. Dotti, Universe 8, 38 (2022), 2201.03640 | arXiv:2201.03640v1; INSPIRE 2008924; Crossref 10.3390/universe8010038 | verified (volume, page and year match in INSPIRE and Crossref) |
| MP2014 | Mbarek, Paranjape, Phys. Rev. D 90, 101502 (2014), 1407.1457 | arXiv:1407.1457v2; INSPIRE 1305085; Crossref 10.1103/PhysRevD.90.101502 | corrected, title, the APS title reads "spacetime" (arXiv: "space-time"); the article is a Rapid Communication, 101502(R) |
| Johnson2020 | Johnson, Paranjape, Savard, Tapia-Arellano, Gen. Relativ. Gravit. 52, 80 (2020), 1910.01774 | arXiv:1910.01774v4; INSPIRE 1757530; Crossref 10.1007/s10714-020-02732-9 | corrected, journal reference added, the entry gave arXiv only |
| Alves2025 | N. de Aguiar Alves, Landulfo, Arderucio Costa, Phys. Rev. D 111, 044027 (2025), 2408.00154 | arXiv:2408.00154v3; INSPIRE 2813778; Crossref 10.1103/PhysRevD.111.044027 | corrected, first initial, the author is Níckolas de Aguiar Alves, not "M." |
| Moinvaziri2026b | N. Moinvaziri, "Stability of a negative-mass boson star," companion paper | paper/2-boson-star/main.tex `\title` | this repository's Paper II, title matches its `\title`, no external source exists |
| Zerilli1970 | F. J. Zerilli, Phys. Rev. Lett. 24, 737 (1970); Phys. Rev. D 2, 2141 (1970) | Crossref 10.1103/PhysRevLett.24.737; Crossref 10.1103/PhysRevD.2.2141 | corrected, the PRD paper has a different title ("Gravitational field of a particle falling in a Schwarzschild geometry analyzed in tensor harmonics"), now given; both volumes, pages and years match |
| ReedSimon2 | Reed, Simon, Methods of Modern Mathematical Physics II (Academic Press, 1975) | Open Library, ISBN 0125850026 (Academic Press, 28 Sep 1975) | verified |
| Weidmann1987 | J. Weidmann, Spectral Theory of Ordinary Differential Operators, LNM 1258 (Springer, 1987) | Crossref 10.1007/BFb0077960; zbMATH ("Lecture Notes in Mathematics, 1258. Berlin etc.: Springer-Verlag. VI, 303 p. (1987)") | verified |
| BirmanSolomyak | Birman, Solomyak, Spectral Theory of Self-Adjoint Operators in Hilbert Space (Reidel, 1987) | Open Library, ISBN 9027721793 (D. Reidel, Dordrecht, 1987) | verified |
| GST1996 | Gesztesy, Simon, Teschl, Amer. J. Math. 118, 571 (1996) | Crossref 10.1353/ajm.1996.0024 | verified (volume, page and year match in Crossref) |
| Giveon2004 | Giveon, Kol, Ori, Sever, JHEP 08, 014 (2004), hep-th/0401209 | arXiv:hep-th/0401209v2; INSPIRE 643424; Crossref 10.1088/1126-6708/2004/08/014; journal volume is the month, 2004-08, article 014 (INSPIRE: JHEP 08 (2004) 014) | verified |

### paper/2-boson-star/main.tex (Paper II)

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| Bondi1957 | as in Paper I | Crossref 10.1103/RevModPhys.29.423 | verified (volume, page and year match in Crossref) |
| Shirokov2026 | N. M. Shirokov, I. Nachevsky, arXiv:2608.24577 | arXiv:2608.24577v2; INSPIRE 3195733 | corrected, authors and title hyphen, as in Paper I |
| GD2006 | as in Paper I | arXiv:gr-qc/0604021v2; INSPIRE 713860; Crossref 10.1088/0264-9381/23/15/021 | verified (volume, page and year match in INSPIRE and Crossref) |
| DG2009 | as in Paper I | arXiv:0809.3615v3; INSPIRE 797174; Crossref 10.1088/0264-9381/26/21/215002 | verified (volume, page and year match in INSPIRE and Crossref) |
| Moinvaziri2026a | N. Moinvaziri, "The instability of negative-mass Schwarzschild is a point-mass effect: exact thresholds for finite stars," companion paper | paper/1-thresholds/main.tex `\title` | this repository's Paper I, title matches its `\title`, no external source exists |
| GW1989 | M. Gleiser, R. Watkins, Nucl. Phys. B 319, 733 (1989) | INSPIRE 266679; Crossref 10.1016/0550-3213(89)90627-5 | verified (volume, page and year match in Crossref) |
| LeePang1989 | T. D. Lee, Y. Pang, Nucl. Phys. B 315, 477 (1989) | INSPIRE 24641; Crossref 10.1016/0550-3213(89)90365-9 | verified (volume, page and year match in Crossref) |
| KMS1991 | Kusmartsev, Mielke, Schunck, Phys. Rev. D 43, 3895 (1991) | INSPIRE 27886; Crossref 10.1103/PhysRevD.43.3895 | verified (volume, page and year match in Crossref) |
| Bonnor1989 | W. B. Bonnor, Gen. Relativ. Gravit. 21, 1143 (1989) | Crossref 10.1007/BF00763458 | verified (volume, page and year match in Crossref) |
| CJM2004 | Cline, Jeon, Moore, Phys. Rev. D 70, 043543 (2004), hep-ph/0311312 | arXiv:hep-ph/0311312v4; INSPIRE 633909; Crossref 10.1103/PhysRevD.70.043543 | verified (volume, page and year match in INSPIRE and Crossref) |
| Kain2021 | B. Kain, "Boson stars and their radial oscillations," Phys. Rev. D 103, 123003 (2021) | INSPIRE 1866961; Crossref 10.1103/PhysRevD.103.123003; arXiv:2106.01740 (INSPIRE record) | corrected, title, "Radial oscillations and stability of multiple-field boson stars" is not the title of PRD 103, 123003 (Kain's "...multiple-fluid compact stars" is PRD 102, 023001) |
| Coleman1985 | S. Coleman, Nucl. Phys. B 262, 263 (1985) | INSPIRE 214529; Crossref 10.1016/0550-3213(85)90286-X | verified (volume, page and year match in Crossref) |
| FLP1987 | Friedberg, Lee, Pang, Phys. Rev. D 35, 3640 (1987) | INSPIRE 20906; Crossref 10.1103/PhysRevD.35.3640 | verified (volume, page and year match in Crossref) |
| Sorkin1982 | R. D. Sorkin, Astrophys. J. 257, 847 (1982) | Crossref 10.1086/160034 | verified (volume, page and year match in Crossref) |
| Alves2025 | N. de Aguiar Alves et al., Phys. Rev. D 111, 044027 (2025), 2408.00154 | arXiv:2408.00154v3; INSPIRE 2813778; Crossref 10.1103/PhysRevD.111.044027 | corrected, first initial, as in Paper I |
| MisnerSharp | Misner, Sharp, Phys. Rev. 136, B571 (1964) | INSPIRE 33605; Crossref 10.1103/PhysRev.136.B571 | verified (volume, page and year match in Crossref) |
| SeidelSuen1990 | Seidel, Suen, Phys. Rev. D 42, 384 (1990) | INSPIRE 27230; Crossref 10.1103/PhysRevD.42.384 | verified (volume, page and year match in Crossref) |
| HawleyChoptuik2000 | Hawley, Choptuik, Phys. Rev. D 62, 104024 (2000) | INSPIRE 530431; Crossref 10.1103/PhysRevD.62.104024 | verified (volume, page and year match in Crossref) |
| Macedo2013 | Macedo, Pani, Cardoso, Crispino, Phys. Rev. D 88, 064046 (2013), 1307.4812 | arXiv:1307.4812v3; INSPIRE 1243407; Crossref 10.1103/PhysRevD.88.064046 | verified (volume, page and year match in INSPIRE and Crossref) |

### paper/3-wormhole-largeD/main.tex (Paper III)

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| MT1988 | Morris, Thorne, Am. J. Phys. 56, 395 (1988) | Crossref 10.1119/1.15620 | verified (volume, page and year match in Crossref) |
| Ellis1973 | H. G. Ellis, J. Math. Phys. 14, 104 (1973) | Crossref 10.1063/1.1666161; an erratum, J. Math. Phys. 15, 520 (1974), exists and is not cited | verified (volume, page and year match in Crossref) |
| Bronnikov1973 | K. A. Bronnikov, Acta Phys. Pol. B 4, 251 (1973) | INSPIRE 87389; INSPIRE: Acta Phys. Polon. B 4 (1973) 251-266, no DOI | verified |
| GGS2009a | González, Guzmán, Sarbach, Class. Quantum Grav. 26, 015010 (2009), 0806.0608 | arXiv:0806.0608v2; INSPIRE 787230; Crossref 10.1088/0264-9381/26/1/015010 | verified (volume, page and year match in INSPIRE and Crossref) |
| SH2002 | Shinkai, Hayward, Phys. Rev. D 66, 044005 (2002), gr-qc/0205041 | arXiv:gr-qc/0205041v2; INSPIRE 586606; Crossref 10.1103/PhysRevD.66.044005 | verified (volume, page and year match in INSPIRE and Crossref) |
| GGS2009b | González, Guzmán, Sarbach, Class. Quantum Grav. 26, 015011 (2009), 0806.1370 | arXiv:0806.1370v2; INSPIRE 787653; Crossref 10.1088/0264-9381/26/1/015011 | verified (volume, page and year match in INSPIRE and Crossref) |
| TS2013 | Torii, Shinkai, Phys. Rev. D 88, 064027 (2013), 1309.2058 | arXiv:1309.2058v1; INSPIRE 1253353; Crossref 10.1103/PhysRevD.88.064027 | verified (volume, page and year match in INSPIRE and Crossref) |
| EST2013 | Emparan, Suzuki, Tanabe, JHEP 06, 009 (2013), 1302.6382 | arXiv:1302.6382v3; INSPIRE 1221232; Crossref 10.1007/JHEP06(2013)009; journal volume is the month, 2013-06, article 009 | verified |
| PSLQ1999 | Ferguson, Bailey, Arno, Math. Comp. 68, 351 (1999) | Crossref 10.1090/S0025-5718-99-00995-3 | verified (volume, page and year match in Crossref) |
| GMG2007 | Gundlach, Martín-García, Living Rev. Relativ. 10, 5 (2007), 0711.4620 | arXiv:0711.4620v1; INSPIRE 768973; Crossref 10.12942/lrr-2007-5 | verified (volume, page and year match in INSPIRE and Crossref) |
| ST2017 | Shinkai, Torii, Phys. Rev. D 96, 044009 (2017), 1706.02070 | arXiv:1706.02070v3; INSPIRE 1603086; Crossref 10.1103/PhysRevD.96.044009 | verified (volume, page and year match in INSPIRE and Crossref) |
| XCY2025 | Xu, Chew, Yeom, JCAP 08 (2025) 012, 2503.07610 | arXiv:2503.07610v1; INSPIRE 2898811; Crossref 10.1088/1475-7516/2025/08/012; JCAP issue 08, article 012 (INSPIRE and DOI 10.1088/1475-7516/2025/08/012) | verified (volume, page and year match in INSPIRE) |

### notes/novelty/RESULTS.md

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| GHI | Gibbons, Hartnoll, Ishibashi, PTP 113, 963 (2005), hep-th/0409307 | arXiv:hep-th/0409307v1; INSPIRE 660779; Crossref 10.1143/PTP.113.963 | corrected, "(2005)" added; the label year 2004 is the arXiv year |
| GD06 | Gleiser, Dotti, CQG 23, 5063, gr-qc/0604021 | arXiv:gr-qc/0604021v2; INSPIRE 713860; Crossref 10.1088/0264-9381/23/15/021 | verified, flagged in "Limits of this search" |
| DG08 | Dotti, Gleiser, CQG 26, 215002 (2009), 0809.3615 | arXiv:0809.3615v3; INSPIRE 797174; Crossref 10.1088/0264-9381/26/21/215002 | verified, flagged in "Limits of this search" |
| CC06 | Cardoso, Cavaglià, PRD 74, 024027, gr-qc/0604101 | arXiv:gr-qc/0604101v2; INSPIRE 715170; Crossref 10.1103/PhysRevD.74.024027 | verified, flagged in "Limits of this search" |
| Dotti2022 | Dotti, Universe 8, 38, 2201.03640 | arXiv:2201.03640v1; INSPIRE 2008924; Crossref 10.3390/universe8010038 | verified (volume, page and year match in INSPIRE and Crossref) |
| Dotti+08 | Dotti, Gleiser, Pullin, Ranea-Sandoval, Vucetich, arXiv:0810.0025 | arXiv:0810.0025v1; INSPIRE 798124; Crossref 10.1142/S0217751X09045030; author list matches; published in Int. J. Mod. Phys. A 24, 1578 (2009), not cited here | verified |
| Giveon04 | Giveon, Kol, Ori, Sever, JHEP 0408:014, hep-th/0401209 | arXiv:hep-th/0401209v2; INSPIRE 643424; Crossref 10.1088/1126-6708/2004/08/014 | verified |
| MP14 | Mbarek, Paranjape, PRD 90, 101502, 1407.1457 | arXiv:1407.1457v2; INSPIRE 1305085; Crossref 10.1103/PhysRevD.90.101502 | verified, flagged in "Limits of this search" |
| Johnson20 | Johnson, Paranjape, Savard, Tapia-Arellano, GRG 52, 80, 1910.01774 | arXiv:1910.01774v4; INSPIRE 1757530; Crossref 10.1007/s10714-020-02732-9 | corrected, page 80 added |
| Novikov18 | Novikov, Bisnovatyi-Kogan, Novikov, PRD 98, 063528 (2018), 1807.06468 | arXiv:1807.06468v1; INSPIRE 1682851; Crossref 10.1103/PhysRevD.98.063528 | verified (volume, page and year match in INSPIRE and Crossref) |
| Alves25 | de Aguiar Alves, Landulfo, Arderucio Costa, PRD 111, 044027 (2025), 2408.00154 | arXiv:2408.00154v3; INSPIRE 2813778; Crossref 10.1103/PhysRevD.111.044027 | verified (volume, page and year match in INSPIRE and Crossref) |
| Paranjape25 | Paranjape, "Puzzles with Negative Mass," Springer (Proc. Gravity, Strings and Fields), pp. 167-173, ISBN 9783031912665 | Crossref 10.1007/978-3-031-91266-5_14; Crossref: CRM Series in Mathematical Physics, pp. 167-173, ISBNs 9783031912658 / 9783031912665 | verified |
| Araneda17 | Araneda, Dotti, PRD 96, 104020 (2017), 1611.03534 | arXiv:1611.03534v3; INSPIRE 1497549; Crossref 10.1103/PhysRevD.96.104020 | verified (volume, page and year match in INSPIRE and Crossref) |
| Kimura17 | Kimura, CQG, 1706.01447 | arXiv:1706.01447v3; INSPIRE 1602611; Crossref 10.1088/1361-6382/aa903f; CQG 34, 235007 (2017), not cited in full | verified |
| Shirokov26 | Shirokov, arXiv:2608.24577 | arXiv:2608.24577v2; INSPIRE 3195733 | verified, the record now lists Shirokov and Nachevsky (v2, 2026-10-26); the note dates from 2026-10-22, when only v1 existed, so it is left unchanged |
| Dzh14 | Dzhunushaliev et al., PRD 90, 124038 (2014), 1409.6978 | arXiv:1409.6978v1; INSPIRE 1318825; Crossref 10.1103/PhysRevD.90.124038 | verified (volume, page and year match in INSPIRE and Crossref) |
| Dzh08 | Dzhunushaliev, Folomeev, Myrzakulov, Singleton, JHEP 0807:094, 0805.3211 | arXiv:0805.3211v3; INSPIRE 786336; Crossref 10.1088/1126-6708/2008/07/094; JHEP 07 (2008) 094 | verified, flagged in "Limits of this search" |
| Dzh16 | Dzhunushaliev et al., PRD 94, 024004 (2016), 1606.07304 | arXiv:1606.07304v2; INSPIRE 1471829; Crossref 10.1103/PhysRevD.94.024004 | verified (volume, page and year match in INSPIRE and Crossref) |
| Bernal26 | Bernal et al., arXiv:2603.04400 | arXiv:2603.04400v1; INSPIRE 3125738; Crossref 10.1103/6qvb-v4q4; published as PRD 113, 104041 (2026), 10.1103/6qvb-v4q4 | verified |
| Jaramillo23 | Jaramillo, Jiménez-Vázquez, Núñez, PRD 107, 064061 (2023), 2302.07889 | arXiv:2302.07889v1; INSPIRE 2633039; Crossref 10.1103/PhysRevD.107.064061 | verified (volume, page and year match in INSPIRE and Crossref) |
| Tamaki11 | Tamaki, Sakai, PRD 81, 124041 (2010), 1105.1498 | arXiv:1105.1498v1; INSPIRE 860321; Crossref 10.1103/PhysRevD.81.124041 | corrected, journal reference added; the label year 2011 is the arXiv posting year, the paper appeared in 2010 |
| Sorkin82 | Sorkin, ApJ 257, 847 (1982) | Crossref 10.1086/160034 | verified (volume, page and year match in Crossref) |
| Ivanov12 | Ivanov, Phys. Usp. 55, 1232-1238 (2012) [UFN 182, 1319], 1204.1269 | arXiv:1204.1269v2; INSPIRE 1104736; Crossref 10.3367/UFNe.0182.201212d.1319; DOI 10.3367/UFNe.0182.201212d.1319 encodes UFN 182, 1319 | verified (volume, page and year match in INSPIRE and Crossref) |
| Forward90 | Forward, J. Propulsion Power 6, 28 (1990) | Crossref 10.2514/3.23219 | verified (volume, page and year match in Crossref) |
| Bondi57 | Bondi, RMP 29, 423 | Crossref 10.1103/RevModPhys.29.423 | verified (volume, page and year match in Crossref) |
| BK09 | Bičák, Kofroň, GRG 41, 153 (2009), 0901.3971 | arXiv:0901.3971v1; INSPIRE 811803; Crossref 10.1007/s10714-008-0662-0 | verified, flagged in "Limits of this search" |
| BS64 | Bonnor, Swaminarayan, Z. Phys. 177, 240 (1964) | Crossref 10.1007/BF01375497; Crossref: pp. 240-256 | verified |
| TL26 | Trivedi, Loeb, arXiv:2605.10976 | arXiv:2605.10976v1; INSPIRE 3154777; Crossref 10.1016/j.dark.2026.102409; published as Phys. Dark Univ. 53, 102409 (2026) | verified |
| Price93 | Price, Am. J. Phys. 61, 216 (1993) | Crossref 10.1119/1.17293 | verified (volume, page and year match in Crossref) |
| PdA14 | Petit, d'Agostini, arXiv:1408.2451 | arXiv:1408.2451v1; arXiv only, not in INSPIRE | verified |
| Cao12 | Cao, IJMPD 21, 1250061 (2012), DOI 10.1142/S0218271812500617 | Crossref 10.1142/S0218271812500617 | verified (volume, page and year match in Crossref) |
| SN26 | Shirokov, Nachevsky, arXiv:2610.13206 | arXiv:2610.13206v2; INSPIRE 3213722; title and authors match | verified |
| CCP26 | Caridi, Corelli, Pani, arXiv:2603.13210 | arXiv:2603.13210v2; INSPIRE 3129166; Crossref 10.1103/3c7y-437m; v2 of 2026-11-06 accepted by PRD, DOI 10.1103/3c7y-437m | verified |
| CLZ24 | Cecchini, Lesourd, Zeidler, IMRN (2024), 2307.05277 | arXiv:2307.05277v2; INSPIRE 2676030; Crossref 10.1093/imrn/rnad315; IMRN 2024(9), 7870-7890 | verified |
| ShiTam04 | Shi, Tam, math/0406559 | arXiv:math/0406559v3; arXiv only, not in INSPIRE; title "Some lower estimates of ADM mass and Brown-York mass" | verified |
| arXiv ids | 2509.05244, 2604.24746, 2508.17585, 2606.23529, 2411.00933, 1406.0841, 2604.00071 | arXiv API: all exist, the topics named in the note match the titles | verified, 2411.00933 is Thong and Melatos (PRD 110, 124010); 1406.0841 is Babichev, Dokuchaev, Eroshenko (Phys. Usp. 56, 1155); 2604.00071 is Shirokov's own wormhole paper, so "cited by Shirokov" there is a self-citation |

### notes/critical_ghost.md

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| SH02 | Shinkai, Hayward, PRD 66, 044005 (2002) | arXiv:gr-qc/0205041v2; INSPIRE 586606; Crossref 10.1103/PhysRevD.66.044005 | verified (volume, page and year match in INSPIRE and Crossref) |
| GGS09a | González, Guzmán, Sarbach, CQG 26, 015010 (2009) | arXiv:0806.0608v2; INSPIRE 787230; Crossref 10.1088/0264-9381/26/1/015010 | verified (volume, page and year match in INSPIRE and Crossref) |
| GGS09b | González, Guzmán, Sarbach, CQG 26, 015011 (2009) | arXiv:0806.1370v2; INSPIRE 787653; Crossref 10.1088/0264-9381/26/1/015011 | verified (volume, page and year match in INSPIRE and Crossref) |
| GGS09c | González, Guzmán, Sarbach, "Instability of charged wormholes supported by a ghost scalar field," PRD 80, 024023 (2009) | arXiv:0906.0420v2; INSPIRE 822042; Crossref 10.1103/PhysRevD.80.024023 | corrected, *(unverified)* removed; title changed to the published one (arXiv: "On the instability of ...") |
| DHNS09 | Doroshkevich, Hansen, Novikov, Shatskiy, IJMPD 18, 1665 (2009) | arXiv:0812.0702v2; INSPIRE 804234; Crossref 10.1142/S0218271809015230 | verified, *(unverified)* removed; arXiv:0811.2558 is marked withdrawn on arXiv and INSPIRE; the arXiv journal-ref field misprints the journal, INSPIRE and the DOI say IJMPD |
| TS13 | Torii, Shinkai, PRD 88, 064027 (2013) | arXiv:1309.2058v1; INSPIRE 1253353; Crossref 10.1103/PhysRevD.88.064027 | verified (volume, page and year match in INSPIRE and Crossref) |
| ST17 | Shinkai, Torii, PRD 96, 044009 (2017) | arXiv:1706.02070v3; INSPIRE 1603086; Crossref 10.1103/PhysRevD.96.044009 | verified, *(unverified)* removed; the volume and article number are filled in |
| BFZ11 | Bronnikov, Fabris, Zhidenko, EPJC 71, 1791 (2011) | arXiv:1109.6576v2; INSPIRE 929989; Crossref 10.1140/epjc/s10052-011-1791-2 | verified, *(unverified)* removed |
| BKZ12 | Bronnikov, Konoplya, Zhidenko, PRD 86, 024028 (2012) | arXiv:1205.2224v3; INSPIRE 1114144; Crossref 10.1103/PhysRevD.86.024028 | verified, *(unverified)* removed |
| CPS20 | Cremona, Pizzocchero, Sarbach, PRD 101, 104061 (2020) | arXiv:1911.13103v2; INSPIRE 1767655; Crossref 10.1103/PhysRevD.101.104061 | verified (volume, page and year match in INSPIRE and Crossref) |
| DFKKR13 | Dzhunushaliev, Folomeev, Kleihaus, Kunz, Radu, PRD 88, 124028 (2013) | arXiv:1309.2448v1; INSPIRE 1253547; Crossref 10.1103/PhysRevD.88.124028 | verified, *(unverified)* replaced by the journal reference |
| DFKK13 | Dzhunushaliev, Folomeev, Kleihaus, Kunz, PRD 87, 104036 (2013) | arXiv:1302.5217v2; INSPIRE 1220529; Crossref 10.1103/PhysRevD.87.104036 | verified (volume, page and year match in INSPIRE and Crossref) |
| Kai26 | B. Kain, "Dynamical wormholes," arXiv:2602.18231 | arXiv:2602.18231v1; INSPIRE 3122144; Crossref 10.1088/1361-6382/ae44e2; published as CQG 43, 045014 (2026) | verified |
| Shi26a | N. M. Shirokov, arXiv:2604.00071 | arXiv:2604.00071v1; INSPIRE 3138462; single author confirmed | verified |
| CCP26 | Caridi, Corelli, Pani, arXiv:2603.13210 | arXiv:2603.13210v2; INSPIRE 3129166; Crossref 10.1103/3c7y-437m; v2 accepted by PRD, DOI 10.1103/3c7y-437m | verified |
| NRM12 | Nakonieczna, Rogatko, Moderski, PRD 86, 044043 (2012) | arXiv:1209.1203v1; INSPIRE 1184618; Crossref 10.1103/PhysRevD.86.044043 | verified (volume, page and year match in INSPIRE and Crossref) |
| HY10 | Hwang, Yeom, CQG 27, 205002 (2010) | arXiv:1002.4246v2; INSPIRE 846339; Crossref 10.1088/0264-9381/27/20/205002 | verified (volume, page and year match in INSPIRE and Crossref) |
| GBH19 | Gundlach, Baumgarte, Hilditch, PRD 100, 104010 (2019) | arXiv:1908.05971v2; INSPIRE 1749935; Crossref 10.1103/PhysRevD.100.104010 | verified (volume, page and year match in INSPIRE and Crossref) |
| Mae09 | H. Maeda, PRD 79, 024030 (2009) | arXiv:0811.2962v2; INSPIRE 802900; Crossref 10.1103/PhysRevD.79.024030 | verified (volume, page and year match in INSPIRE and Crossref) |
| Fro97 | A. V. Frolov, PRD 56, 6433 (1997) | arXiv:gr-qc/9704040v1; INSPIRE 442210; Crossref 10.1103/PhysRevD.56.6433 | verified, *(unverified)* removed |
| HK04 | Hayward, Koyama, PRD 70, 101502 (2004) | arXiv:gr-qc/0406080v1; INSPIRE 652760; Crossref 10.1103/PhysRevD.70.101502 | verified, *(unverified)* replaced by the journal reference |
| Ber26 | Bernal et al., PRD 113, 104041 (2026) | arXiv:2603.04400v1; INSPIRE 3125738; Crossref 10.1103/6qvb-v4q4 | verified, *(unverified)* replaced by the journal reference |
| JJN23 | Jaramillo, Jiménez-Vázquez, Núñez, "Confinement of exotic matter: static solutions," PRD 107, 064061 (2023) | arXiv:2302.07889v1; INSPIRE 2633039; Crossref 10.1103/PhysRevD.107.064061 | corrected, initial "E." added to Jiménez-Vázquez; title changed to the published one (arXiv adds "I."); "as recorded in notes/novelty" dropped |
| DMW23 | Ding, Ma, Fang, Wang, JHEP 04 (2024) 033 | arXiv:2305.19819v1; INSPIRE 2664190; Crossref 10.1007/JHEP04(2024)033; JHEP 04 (2024) 033 | corrected, *(unverified)* removed; the published paper has four authors (T.-F. Fang added), arXiv v1 has three |
| CDK24 | Clough, Dietrich, Khan, doi:10.33232/001c.121868 (Open J. Astrophys., 2024) | arXiv:2406.02466v2; INSPIRE 2794423; Crossref 10.33232/001c.121868; Open J. Astrophys. 7 (2024) | verified |
| GHMG25 | Gundlach, Hilditch, Martín-García, arXiv:2507.07636 | arXiv:2507.07636v1; INSPIRE 2944628; authors match | verified |
| LRR07 | Gundlach, Martín-García, LRR 10, 5 (2007) | arXiv:0711.4620v1; INSPIRE 768973; Crossref 10.12942/lrr-2007-5 | verified (volume, page and year match in INSPIRE and Crossref) |
| GMG96 | Gundlach, Martín-García, PRD 54, 7353 (1996) | arXiv:gr-qc/9606072v2; INSPIRE 420116; Crossref 10.1103/PhysRevD.54.7353 | verified (volume, page and year match in INSPIRE and Crossref) |
| MGM26 | Mittal, Gundlach, Martel, arXiv:2609.17207 | arXiv:2609.17207v1; INSPIRE 3204115; authors match, no journal reference yet | verified |
| LC96 | Liebling, Choptuik, PRL 77, 1424 (1996) | arXiv:gr-qc/9606057v1; INSPIRE 419789; Crossref 10.1103/PhysRevLett.77.1424 | verified (volume, page and year match in INSPIRE and Crossref) |
| HE97 | Hirschmann, Eardley, PRD 56, 4696 (1997) | arXiv:gr-qc/9511052v1; INSPIRE 402477; Crossref 10.1103/PhysRevD.56.4696 | verified (volume, page and year match in INSPIRE and Crossref) |
| HHS96 | Hamadé, Horne, Stewart, CQG 13, 2241 (1996) | arXiv:gr-qc/9511024v2; INSPIRE 401928; Crossref 10.1088/0264-9381/13/8/017 | verified (volume, page and year match in INSPIRE and Crossref) |
| EHH95 | Eardley, Hirschmann, Horne, PRD 52, R5397 (1995) | arXiv:gr-qc/9505041v1; INSPIRE 395318; Crossref 10.1103/PhysRevD.52.R5397 | verified (volume, page and year match in INSPIRE and Crossref) |
| AH20b | Antonelli, Hatefi, JHEP 03 (2020) 180 | arXiv:1912.06103v3; INSPIRE 1770477; Crossref 10.1007/JHEP03(2020)180; JHEP 03 (2020) 180 | verified |
| HK21 | Hatefi, Kuntz, EPJC 81 (2021) | arXiv:2010.11603v4; INSPIRE 1824435; Crossref 10.1140/epjc/s10052-020-08788-w | verified (volume, page and year match in INSPIRE and Crossref) |
| LTHA02 | Lechner, Thornburg, Husa, Aichelburg, PRD 65, 081501 (2002) | arXiv:gr-qc/0112008v1; INSPIRE 567911; Crossref 10.1103/PhysRevD.65.081501 | verified (volume, page and year match in INSPIRE and Crossref) |
| CCB96 | Choptuik, Chmaj, Bizoń, PRL 77, 424 (1996) | arXiv:gr-qc/9603051v1; INSPIRE 417159; Crossref 10.1103/PhysRevLett.77.424 | verified (volume, page and year match in INSPIRE and Crossref) |
| CHM99 | Choptuik, Hirschmann, Marsa, PRD 60, 124011 (1999) | arXiv:gr-qc/9903081v1; INSPIRE 497094; Crossref 10.1103/PhysRevD.60.124011 | verified (volume, page and year match in INSPIRE and Crossref) |
| BCG97 | Brady, Chambers, Gonçalves, PRD 56, R6057 (1997) | arXiv:gr-qc/9709014v1; INSPIRE 448088; Crossref 10.1103/PhysRevD.56.R6057 | verified (volume, page and year match in INSPIRE and Crossref) |

### notes/bh_to_wormhole.md

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| H94 | Hayward, PRD 49, 6467 (1994), gr-qc/9303006 | arXiv:gr-qc/9303006v3; INSPIRE 353003; Crossref 10.1103/PhysRevD.49.6467; the published title is "General laws of black-hole dynamics"; the INSPIRE record carries the v1 title "Marginal surfaces and apparent horizons" | verified (volume, page and year match in Crossref) |
| H99 | Hayward, IJMPD 8, 373 (1999), gr-qc/9805019 | arXiv:gr-qc/9805019v3; INSPIRE 470204; Crossref 10.1142/S0218271899000286 | verified, "(journal from citations)" removed, the journal is confirmed directly |
| H02 | Hayward, gr-qc/0203051 (conference contribution; journal not found) | arXiv:gr-qc/0203051v1; INSPIRE 584097 | verified, INSPIRE lists it as a conference paper (cnum C02-01-09) with no journal reference, so "journal not found" stands |
| H09 | Hayward, PRD 79, 124001 (2009), 0903.5438 | arXiv:0903.5438v1; INSPIRE 816798; Crossref 10.1103/PhysRevD.79.124001 | verified (volume, page and year match in INSPIRE and Crossref) |
| HK04 | Hayward, Koyama, PRD 70, 101502 (2004) | arXiv:gr-qc/0406080v1; INSPIRE 652760; Crossref 10.1103/PhysRevD.70.101502 | verified, "(journal from [H09]'s reference list)" removed, confirmed directly |
| KH04 | Koyama, Hayward, PRD 70, 084001 (2004) | arXiv:gr-qc/0406113v2; INSPIRE 653241; Crossref 10.1103/PhysRevD.70.084001 | verified (volume, page and year match in INSPIRE and Crossref) |
| HKL01 | Hayward, Kim, Lee, PRD 65, 064003 (2002), gr-qc/0110080 (2D dilaton gravity) | arXiv:gr-qc/0110080v1; INSPIRE 564565; Crossref 10.1103/PhysRevD.65.064003 | corrected, *(unverified)* replaced by the journal reference; the abstract confirms the 2D (CGHS) dilaton gravity description |
| GMG96 | Gundlach, Martín-García, PRD 54, 7353 (1996) | arXiv:gr-qc/9606072v2; INSPIRE 420116; Crossref 10.1103/PhysRevD.54.7353 | verified (volume, page and year match in INSPIRE and Crossref) |
| SMMV19 | Simpson, Martín-Moruno, Visser, CQG 36, 145007 (2019) | arXiv:1902.04232v3; INSPIRE 1719970; Crossref 10.1088/1361-6382/ab28a5 | verified (volume, page and year match in INSPIRE and Crossref) |
| CK21 | Chakrabarti, Kar, PRD 104, 024071 (2021) | arXiv:2106.14761v1; INSPIRE 1870562; Crossref 10.1103/PhysRevD.104.024071 | verified (volume, page and year match in INSPIRE and Crossref) |
| Ter22 | D. R. Terno, PRD 106, 044035 (2022), 2203.03770 | arXiv:2203.03770v2; INSPIRE 2048052; Crossref 10.1103/PhysRevD.106.044035 | corrected, "journal not found" replaced by the journal reference |
| NM08 | Nozawa, Maeda, CQG 25, 055009 (2008), 0710.2709 | arXiv:0710.2709v3; INSPIRE 764151; Crossref 10.1088/0264-9381/25/5/055009 | verified (volume, page and year match in INSPIRE and Crossref) |
| MGM26 | Mittal, Gundlach, Martel, arXiv:2609.17207 | arXiv:2609.17207v1; INSPIRE 3204115; authors match | verified |

### notes/wormhole_largeD/one_mode.md

| Key | Entry as cited now | Source used | Result |
|---|---|---|---|
| GGS09a | González, Guzmán, Sarbach, CQG 26, 015010 (2009) | arXiv:0806.0608v2; INSPIRE 787230; Crossref 10.1088/0264-9381/26/1/015010 | verified (volume, page and year match in INSPIRE and Crossref) |
| CPS20 | Cremona, Pizzocchero, Sarbach, PRD 101, 104061 (2020) | arXiv:1911.13103v2; INSPIRE 1767655; Crossref 10.1103/PhysRevD.101.104061 | verified (volume, page and year match in INSPIRE and Crossref) |
| TS13 | Torii, Shinkai, PRD 88, 064027 (2013) | arXiv:1309.2058v1; INSPIRE 1253353; Crossref 10.1103/PhysRevD.88.064027 | verified (volume, page and year match in INSPIRE and Crossref) |
| Wei87 | Weidmann, LNM 1258 (Springer, Berlin, 1987) | Crossref 10.1007/BFb0077960; zbMATH ("Lecture Notes in Mathematics, 1258. Berlin etc.: Springer-Verlag. VI, 303 p. (1987)") | verified, *(unverified)* removed |
| BS91 | Berezin, Shubin, The Schrödinger Equation (Springer/Kluwer, Dordrecht, 1991) | Crossref 10.1007/978-94-011-3154-4; Open Library (Kluwer Academic Publishers, 1991) | verified, *(unverified)* removed; the publisher is Kluwer, now distributed by Springer Netherlands (Crossref) |
| Sim05 | Simon, "Sturm oscillation and comparison theorems," in Sturm-Liouville Theory: Past and Present (Birkhäuser, Basel, 2005), pp. 29-43 | Crossref 10.1007/3-7643-7359-8_2; Crossref 10.1007/3-7643-7359-8 (editors Amrein, Hinz, Pearson; Birkhäuser Basel, 2005) | verified, *(unverified)* removed; the page range is added |

# CONTEXT: everything an agent needs to pick up this project

This file explains the whole project from start to finish: what it is, how
it was built, what was decided and why, what was found, and how to work on
it. Read it first, then `PLAN.md` for the current to-do list.

Last updated: 2026-11-24. State: **complete**. The analysis, notebook and
both reports are committed and pushed to `origin/main`
(https://github.com/nanfaweb/Luciferase-Convergent-Evolution).

---

## 1. What the project is

- **Type.** An undergraduate bioinformatics course project (FAST-NU Lahore,
  2026, Weeks 1–15 syllabus). It is a course deliverable, not a publication.
- **Title.** *Molecular Signatures of Convergent Evolution: Comparative
  Sequence Analysis of Independently Evolved Luciferases.*
- **Question.** Bioluminescence evolved independently many times, and each
  time a different enzyme (a luciferase) was recruited. Do these
  independently evolved luciferases share any **protein-sequence** or
  **physicochemical** features?
- **Level of comparison.** Protein level, not DNA. The user originally framed
  it as "shared DNA"; we settled on proteins, because unrelated genes share
  nothing at DNA level and proteins are where function lives.
- **Hypotheses.**
  - **H1:** independent families share no sequence, motifs or peptides.
  - **H2:** luciferases shift in the same direction from their own
    non-luminous relatives (controls) across families.
  - **H0 for H2:** properties follow ancestry, and the shift directions are
    random.
- **Scope contract.** The original proposal is in
  `docs/Bioluminescence Proposal.pdf` and `.docx`. Background and the
  original brief are in `docs/bioinformatics_project_handoff.md`.

## 2. Answer (the result in one paragraph)

Independently evolved luciferases converged in **function, not sequence**.
Across 1,259 pairs of luciferases from unrelated families, median identity
is about 20% and the median Smith–Waterman z-score is 0, which is chance
level. BLAST finds no hits between independent families, and no PROSITE
motif or tryptic peptide is shared between them. The only cross-family
similarity is firefly with click beetle (about 49% identity), which share an
ancestral acyl-CoA synthetase family.

Physicochemical profiles sort by family. Clustering recovers the families
with an adjusted Rand index (ARI) of 0.80, and Kruskal–Wallis p < 1e-5.

In all 5 families that have controls, luciferases are **less hydrophobic
(lower GRAVY) and more aromatic** than their controls. The sign test gives
p = 0.0625 for each, the smallest value possible with 5 families. This is a
consistent hint, **not a significant finding**. Conclusions must stay
cautious: never claim that a property causes light production.

## 3. Dataset

`data/seeds.csv` is the **source of truth**. It has one row per protein,
with these columns:

| Column | Meaning |
|---|---|
| `accession` | Protein ID |
| `source` | `uniprot` or `ncbi` |
| `family` | Luciferase family |
| `role` | `luciferase` or `control` |
| `label`, `organism` | Display name and species |
| `allow_partial` | Keep a near-complete sequence (flagged) |
| `strip_suffix` | Text to remove from the record |
| `reason` | Why this protein was included |

`src/download.py` turns it into `data/curated/`, which is the data to
analyse. `data/raw/` is gitignored; never analyse it.

67 proteins: 55 luciferases in 8 families (9 independent origins), plus 12
controls.

| Family | Luc. | Controls | Notes |
|---|---|---|---|
| firefly | 8 | 2: *Drosophila* CG6178 (Q9VCC6), *Tenebrio* luciferase-like (Q1ET68) | ACS family |
| click_beetle | 4 | 1: AbLL A9ZPM4 (non-luminous click beetle) | ACS family, independent origin from fireflies (He et al. 2024). AbLL was moved here from firefly. |
| renilla_type | 3 | 3: haloalkane dehalogenases (LinB, DhaA, DafA) | 2 *Renilla* + 1 brittle star (AfLuc). These are two independent recruitments of the HLD family (Lau et al. 2025). |
| bacterial_luxa | 10 | 3: SsuD (P80645, Q9HYG2, P40402) | Luciferase-like monooxygenase family |
| fungal_luz | 9 | 3: Luz-domain proteins (A0A0D7AKP7, A0A6A4GH72, A0A550C6L9) | |
| copepod | 10 | – | No InterPro entry, which is expected |
| dinoflagellate | 8 | – | Three-domain enzyme |
| ostracod | 3 | – | Cypridinid |

**Selection rules.**
- There must be functional evidence that the protein is a luciferase.
- Sequences must be complete. Six near-complete clones are kept and flagged
  with `allow_partial`.
- No two sequences within a family may be ≥97% identical.
- Reviewed entries are preferred.
- Records are versioned and hashed. `python -m src.verify_data --live`
  re-checks them against live UniProt.

**Sources.** 62 proteins come from UniProt and 5 from NCBI. Those 5 are the
only proteins without an AlphaFold model. More detail is in
`data/README.md`.

**How the dataset grew** (see the git log):

| Commit | Dataset |
|---|---|
| `e817a62` | 21 luciferases + 6 controls |
| `c895ef7` | 30 + 7, with NCBI sources added |
| `c0aa890` | Copepod, dinoflagellate and ostracod families added (51 + 7 = 58) |
| `eb9b624` | Click beetles and 5 more controls added (55 + 12 = 67) |

## 4. Pipeline: modules, outputs, experiments

Run everything from the project root. **The order matters.**

| # | Command | Course week | What it does → outputs |
|---|---|---|---|
| – | `python -m src.download` then `python -m src.verify_data --live` | W3 | Fetch, check and curate → `data/curated/*.fasta`, `manifest.csv`, `validation_summary.json` |
| E1 | `python -m src.features` | W14 | ProtParam properties: length, MW, pI, GRAVY, aromaticity, `aa_X_pct` → `results/features.csv`, `features_by_family.csv`. Exports `RESULTS_DIR`. |
| E2, E3 | `python -m src.align` (about 3 min) | W6–W8 | All 2,211 pairs: NW identity, SW score with a z-score against 20 shuffles. Uses BLOSUM62, gap −10/−0.5, end gaps 0 for global. Also the hand-written NW DP (`needleman_wunsch_dp`) on HEAGAWGHEE/PAWHEAE, gap −8; its score of −8 matches Biopython. → `pairwise.csv`, `identity_matrix.csv`, `similarity_summary.csv`, `example_alignments.txt`, `dp_example.txt`, `figures/identity_heatmap.png` |
| E4 | `python -m src.dotplot` | W9 | Window 10, threshold 23 → `figures/dotplots.png` |
| E5 | `python -m src.blast` | W12 | Local BLAST+ all-vs-all, E < 1e-3 → `blast_all_vs_all.tsv`, `blast_summary.csv` |
| E6 | `python -m src.phylo` (needs align) | W10 | Within-family NJ trees on Poisson distance, rooted on the control least identical to the luciferases (midpoint if there is no control) → `results/trees/*.nwk`, `figures/trees.png` |
| E9 | `python -m src.analysis` (needs features) | W11 | PCA, Ward clustering vs family (ARI), Kruskal–Wallis, luciferase-vs-own-control sign test → `stats.csv`, `control_comparison.csv`, `figures/pca_features.png`, `feature_boxplots.png` |
| E7 | `python -m src.motifs` | W2, W3, W12 | PROSITE patterns from UniProt cross-references, converted to regex and scanned against every protein. Also a GO summary. → `prosite_patterns.csv`, `motif_hits.csv`, `go_by_family.csv` |
| E8 | `python -m src.digest` | W5, W14 | Trypsin cut `(?<=[KR])(?!P)`, peptides of 6–30 aa, monoisotopic mass → `tryptic_peptides.csv`, `tryptic_summary.csv` |
| E10 | `python -m src.structure` | W15 | AlphaFold DB API, mean pLDDT (`globalMetricValue`) → `alphafold_confidence.csv` |
| – | `python -m src.report_pdf` (last) | – | Builds 3 summary figures (`similarity_by_category`, `control_differences`, `alphafold_plddt`) and `docs/FINAL_REPORT.pdf`. Numbers are read live from `results/`; references come from `docs/FINAL_REPORT.md`. |
| – | `python -m nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=1200 notebooks/01_analysis.ipynb` | – | Re-runs everything; executes with 0 errors |

**Shared helpers:**
- `src/download.py`: `FAMILIES`, `CURATED_DIR`.
- `src/align.py`: `load_proteins()`, which sorts by family, luciferases
  first; `GLOBAL` and `LOCAL` aligners.
- `src/verify_data.py`: `parse_fasta`.
- `src/features.py`: `RESULTS_DIR`.

**Topics deliberately not used**, as stated in the report:
- W1 (origin of replication), W4 (genome assembly), W9 (sorting by
  reversals), W11 (transcriptomics) and W13 (HIV): there is no DNA or
  genome data here.
- HMMs: used only through Pfam/InterPro domains during curation. No HMM was
  trained, because unrelated luciferases cannot be aligned into one profile.

## 5. Key numbers (67-protein run)

**Pairwise alignment**

| Pair type | Pairs | Median identity | Median SW z |
|---|---|---|---|
| Same family | 194 | 68% | 136 |
| Luciferase vs own control | 86 | 36% | 47 |
| Firefly vs click beetle | 32 | 49% | 126 |
| Other cross-family | 1,259 | 20% | 0.0 (36 pairs above z = 3, max 6.9, none backed by BLAST) |

**BLAST, motifs and peptides**
- **BLAST:** 630 hits. The only cross-family hits are firefly ↔ click
  beetle.
- **Motifs:** 2 PROSITE patterns. AMP_BINDING occurs in firefly and
  click-beetle luciferases **and** their controls. The
  BACTERIAL_LUCIFERASE signature occurs in all LuxA proteins and in no SsuD
  control.
- **Peptides:** 1,806 tryptic peptides. Only 2 occur in more than one
  family (firefly and click beetle).

**Properties and structure**
- **ARI:** 0.80.
- **Sign tests (5 families):**

  | Property | Result | p |
  |---|---|---|
  | GRAVY | lower in 5/5 | 0.0625 |
  | Aromaticity | higher in 5/5 | 0.0625 |
  | Cys % | higher in 4/5 | 0.375 |
  | pI | higher in 2/5 | 1.0 |

- **AlphaFold:** 62 of 67 proteins have models. Median pLDDT is about 90 or
  higher for the beetle, bacterial, fungal and *Renilla* families, 86 for
  dinoflagellates and 56 for copepods.
- **Trees:** the brittle-star luciferase groups with its own DafA control,
  not with *Renilla*.
- **Dot plots:** 479 dots for same family, 334 for luciferase vs control
  and 6 for firefly vs LuxA.

## 6. Design decisions (and why)

1. **Matched non-luminous controls.** Comparing a luciferase with a relative
   from its own family cancels out family traits. Only these 5 families have
   plausible controls.
2. **The family is the statistical unit.** Proteins within a family are not
   independent, so the sign test runs over families.
3. **No mixed-family tree or MSA.** Unrelated proteins cannot be
   meaningfully aligned together.
4. **Composition features, not length or weight,** across families. Size
   reflects domain architecture.
5. **Shuffle z-scores and E-values** give every similarity claim a
   chance baseline.
6. **Signal peptides are kept,** which is a known limitation.

## 7. Documents

| File | Purpose |
|---|---|
| `PLAN.md` | Status, Next, Done, Decisions; the original detailed plan below them |
| `CLAUDE.md` | Stack and the run commands |
| `README.md` | Overview and course-concept table |
| `data/README.md` | Dataset details and limitations |
| `results/README.md` | What each output file is |
| `docs/PROJECT_REPORT.md` | Plain-language, phase-by-phase log for students and instructors. Phase numbers are from the 58-protein run, followed by a 2026-11-23 update section. |
| `docs/FINAL_REPORT.md` | Formal report, in the structure the user requires: problem formulation, existing solutions, proposed methodology/framework, implementation/experiments, results, conclusion, references |
| `docs/FINAL_REPORT.pdf` | The same report as a 16-page PDF: 8 figures, 9 tables, one subsection per course topic (4.3.x) |

**References.** There are 20. All were verified on 2026-11-23 against
Crossref, and each has a DOI. If you change them, edit
`docs/FINAL_REPORT.md` and re-run `src.report_pdf`.

## 8. Environment gotchas

- **OS and shells:** Windows 11. Use PowerShell or Git Bash, with Python
  3.13 and Biopython 1.86.
- **BLAST+ and paths:** BLAST+ 2.17.0 is in `~/tools/ncbi-blast-2.17.0+/bin`,
  on the PATH. It **fails on paths with spaces**, and this folder is
  "Bio Proj", so `src/blast.py` works in `tempfile.TemporaryDirectory()`.
- **Notebook:** `jupyter` is not on PATH. Use `python -m nbconvert`.
- **Amino-acid percentages:** Biopython's `amino_acids_percent` returns
  percentages, which is why the columns are named `aa_X_pct`.
- **UniProt flags:** a "Precursor" flag is **not** a fragment. Only
  "Fragment" is.
- **CSV diffs on re-runs:** the CSVs are written with `lineterminator="\n"`.
  Re-runs can produce formatting-only diffs.
- **PDF tools:** `pymupdf` is installed for rendering PDF pages to check
  them; `pypdf` is not.
- **Git:** pushes occasionally return a 500 error; retrying works.

## 9. Working rules (from the user and the workspace)

**From the workspace `CLAUDE.md`**
- Read `PLAN.md` first.
- Stay in scope: new ideas go into Next.
- Run code before calling something done.
- Update `PLAN.md` and commit with a plain message when a piece is done.
- Never push unless asked. In this project the user has asked for a commit
  and push after each milestone.
- Scratch work goes in `tmp/` (gitignored).
- Never read `.env`.
- Never delete a project; archive it.

**From the user**
- Build only the stage the user names. They delegate design choices, not
  scope.
- Keep `docs/PROJECT_REPORT.md` updated in simple language after each phase.
- Keep the user informed of progress.
- The user asked to be told, and for work to pause, when weekly usage
  reaches 75%. It was 67% on 2026-11-23.
- Commit messages end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## 10. Open items (optional, not requested)

- Compare AlphaFold 3D structures, for example the active-site pocket, to
  look for structural convergence.
- Test the GRAVY and aromaticity hint on more families with controls.
- Analyse mature proteins, with signal peptides removed.
- Rename the folder to `bio-proj`, which the workspace naming rule asks for.
  This is left to the user.

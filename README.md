# Convergent Evolution of Luciferases: Comparative Sequence Analysis

Bioluminescence evolved independently many times. Unrelated protein families were each recruited to produce light. This bioinformatics course project compares **62 luminous proteins from 10 families (11 independent origins)** with **19 non-luminous relatives**, to ask whether independently evolved light-producing proteins share sequence-level or physicochemical features.

## Key findings

- **No shared sequence (H1 supported).** Independent families share no detectable sequence similarity, motifs or tryptic peptides. The only exception is firefly and click beetle, which share an acyl-CoA synthetase ancestor.
- **Composition follows ancestry.** Clustering recovers the families with an adjusted Rand index of 0.82.
- **Pre-registered test (novelty).** A prediction ("luminous proteins are less hydrophobic and more aromatic than their own relatives") was written into git before two new families were added.
  - The **replication failed** for hydrophobicity: the Oplophorus shrimp luciferase is more hydrophobic than its relatives.
  - **Aromaticity holds** across 7 families: family-stratified permutation p = 0.010 (Holm 0.021).
- **The aromatic excess is on the surface** (exploratory), using AlphaFold models: p = 0.0001 for exposed residues, p = 0.70 for buried ones.
- **Shared substrate does not produce shared composition.** Four unrelated coelenterazine-using families are not closer in composition than random families (p = 0.995).

These are associations, not causes. See [`docs/FINAL_REPORT.md`](docs/FINAL_REPORT.md) for the full write-up and limitations.

## Read the report

- [`docs/FINAL_REPORT.pdf`](docs/FINAL_REPORT.pdf): full report with figures and tables (23 pages).
- [`docs/FINAL_REPORT.md`](docs/FINAL_REPORT.md): the same report in Markdown, with 24 references checked against Crossref.
- [`docs/PROJECT_REPORT.md`](docs/PROJECT_REPORT.md): plain-language, phase-by-phase log.
- [`CONTEXT.md`](CONTEXT.md): complete hand-over for another agent or contributor.

## Methods

- **Data:** UniProt and NCBI Protein, curated from `data/seeds.csv` with recorded selection rules.
- **Similarity:** Needleman–Wunsch and Smith–Waterman alignment (BLOSUM62), with shuffle z-scores; dot plots; local BLAST+.
- **Phylogeny:** within-family neighbour-joining trees only (no forced cross-family tree).
- **Properties:** ProtParam (GRAVY, aromaticity, pI, amino-acid composition), PCA, Ward clustering, Kruskal–Wallis.
- **Motifs and peptides:** PROSITE multi-pattern scan, in-silico trypsin digest, GO summary.
- **Structure:** AlphaFold confidence (pLDDT) and Cα burial analysis.
- **Novelty tests:** pre-registered sign test and within-family label permutation (10,000 permutations, seed 1), robustness checks, and a calibration on random features.

## Setup

Requires Python 3.13.

```powershell
python -m pip install -r requirements.txt
python -m src.download
python -m src.verify_data --live
python -m src.features
```

Then run the remaining modules in the order given in [`CLAUDE.md`](CLAUDE.md), or execute the notebook `notebooks/01_analysis.ipynb`. To rebuild the PDF after the analyses: `python -m src.report_pdf`.

BLAST+ 2.17.0 must be installed and on the PATH. `src/blast.py` runs it in a temporary folder, because BLAST cannot handle spaces in paths.

## Project layout

```
data/seeds.csv     source of truth for the dataset (one row per protein, with reasons)
data/raw/          downloaded FASTA (gitignored)
data/curated/      cleaned FASTA and manifest
data/structures/   cached AlphaFold models (gitignored)
src/               pipeline modules
notebooks/         end-to-end analysis notebook
results/           tables and figures
docs/              final report (MD and PDF), project log, handoff notes
```

## Course concepts used

| Week | Topic | Where it appears |
|---|---|---|
| W2, W12 | Motif search, multiple pattern matching | `src/motifs.py`, `results/motif_hits.csv` |
| W3 | Biological databases, GO | UniProt downloads, `data/curated/manifest.csv`, `results/go_by_family.csv` |
| W5, W14 | Mass spectrometry, peptides | `src/digest.py`, `results/tryptic_summary.csv` |
| W6 | Similarity and homology | `src/align.py`, `results/figures/identity_heatmap.png` |
| W7–W8 | Global/local alignment, dynamic programming | `src/align.py`, `results/example_alignments.txt`, `results/dp_example.txt` |
| W9 | Dot plots | `src/dotplot.py`, `results/figures/dotplots.png` |
| W10 | Phylogeny | `src/phylo.py`, `results/figures/trees.png` |
| W11 | Clustering and statistics | `src/analysis.py`, `results/figures/pca_features.png` |
| W11 + novelty | Pre-registered, family-stratified convergence tests | `src/convergence.py`, `results/convergence_tests.csv` |
| W12 | BLAST | `src/blast.py`, `results/blast_summary.csv` |
| W14 | Proteomics | `src/features.py`, `results/features.csv` |
| W15 | Deep learning (AlphaFold) | `src/structure.py`, `src/pockets.py`, `results/aromatic_burial.csv` |

W1, W4, W9 (sorting by reversals), W11 (transcriptomics) and W13 (HIV) are not applicable to protein comparison. The report explains why. W13 (HMMs) is used only through Pfam/InterPro domains during curation.

## Stack

Python 3.13, Biopython 1.86, pandas, NumPy, SciPy, Matplotlib, scikit-learn, reportlab, BLAST+ 2.17.0.

## Licence and data

Sequences come from UniProt and NCBI and remain the property of their sources. Code and documentation are for coursework.

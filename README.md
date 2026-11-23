# Comparative Sequence Analysis of Independently Evolved Luciferases

Course project in bioinformatics: compare luciferase proteins from fireflies, *Renilla*-type, luminous bacteria, fungi, copepods, dinoflagellates and ostracods to ask whether independently evolved light-producing proteins share common sequence-level or physicochemical features.

## Question

Do luciferases from different evolutionary origins remain distinct at the sequence level while showing any shared molecular properties?

## Methods

- Public protein sequences from UniProt / NCBI (FASTA)
- Pairwise similarity, global/local alignment, DotPlots, and BLASTP
- Within-family phylogeny only (no forced cross-family tree)
- Physicochemical features (length, MW, GRAVY, aromaticity, pI) and feature clustering

## Project layout

```
data/raw/          downloaded FASTA
data/curated/      cleaned FASTA + manifest
notebooks/         analysis walkthrough
src/               pipeline modules
results/           tables and figures
```

## Setup

```bash
python -m pip install -r requirements.txt
python -m src.download
python -m src.verify_data --live
python -m src.features
```

`python -m src.download` fetches every accession in `data/seeds.csv` from UniProt (or NCBI Protein, where UniProt lacks it) and checks each record. It then writes `data/raw/` and `data/curated/`. `python -m src.verify_data` re-checks the files offline; with `--live` it also checks them against the current UniProt records. The dataset (51 luciferases from 7 families, 7 non-luminous controls), the selection rules and the limitations are documented in [`data/README.md`](data/README.md).

For a plain-language account of methods and results so far, read [`docs/PROJECT_REPORT.md`](docs/PROJECT_REPORT.md).

See [`PLAN.md`](PLAN.md) for scope, syllabus mapping, and execution phases. Background context is in [`docs/bioinformatics_project_handoff.md`](docs/bioinformatics_project_handoff.md).

## Course concepts used

| Week | Topic | Where it appears |
|---|---|---|
| W3 | Biological databases, GO | UniProt downloads, `manifest.csv`, InterPro/GO notes |
| W6 | Similarity / homology | `src/align.py`, `results/figures/identity_heatmap.png` |
| W7–W8 | Global/local alignment, DP | `src/align.py`, `results/example_alignments.txt` |
| W9 | DotPlots | `src/dotplot.py`, `results/figures/dotplots.png` |
| W10 | Phylogeny | `src/phylo.py`, `results/figures/trees.png` |
| W11 | Clustering | `src/analysis.py`, `results/figures/pca_features.png` |
| W12 | BLAST | `src/blast.py`, `results/blast_summary.csv` |
| W14 | Proteomics | `src/features.py`, `results/features.csv` |

W1–W2 (intro/motifs via InterPro), W13 (HMMs as Pfam/InterPro profiles), and W15 (deep learning) are mentioned in presentation notes only.

## Stack

Python, Biopython, Pandas, NumPy, Matplotlib, SciPy, scikit-learn, BLAST+

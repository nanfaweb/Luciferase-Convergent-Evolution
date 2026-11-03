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
```

`python -m src.download` fetches every accession in `data/seeds.csv` from UniProt (or NCBI Protein, where UniProt lacks it) and checks each record. It then writes `data/raw/` and `data/curated/`. `python -m src.verify_data` re-checks the files offline; with `--live` it also checks them against the current UniProt records. The dataset (51 luciferases from 7 families, 7 non-luminous controls), the selection rules and the limitations are documented in [`data/README.md`](data/README.md).

See [`PROJECT_PLAN.md`](PROJECT_PLAN.md) for scope, syllabus mapping, and execution phases. Background context is in [`bioinformatics_project_handoff.md`](bioinformatics_project_handoff.md).

## Course concepts used

| Week | Topic | Where it appears |
|---|---|---|
| W3 | Biological databases, GO | UniProt downloads, `manifest.csv`, InterPro/GO notes |
| W6 | Similarity / homology | Pairwise identity matrix (Phase 2) |
| W7–W8 | Global/local alignment, DP | Needleman–Wunsch / Smith–Waterman (Phase 2) |
| W9 | DotPlots | Within- vs across-family DotPlots (Phase 2) |
| W10 | Phylogeny | Within-family trees only (Phase 2) |
| W11 | Clustering | Physicochemical feature clustering (Phase 2) |
| W12 | BLAST | BLASTP against curated set (Phase 2) |
| W14 | Proteomics | ProtParam features (Phase 2) |

W1–W2 (intro/motifs via InterPro), W13 (HMMs as Pfam/InterPro profiles), and W15 (deep learning) are mentioned in presentation notes only.

## Stack

Python, Biopython, Pandas, NumPy, Matplotlib, SciPy, Seaborn

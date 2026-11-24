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

Then run the analysis modules listed in `CLAUDE.md`, or the notebook `notebooks/01_analysis.ipynb`. The full write-up is [`docs/FINAL_REPORT.md`](docs/FINAL_REPORT.md), also as a PDF with figures and tables: [`docs/FINAL_REPORT.pdf`](docs/FINAL_REPORT.pdf) (rebuild with `python -m src.report_pdf`).

`python -m src.download` fetches every accession in `data/seeds.csv` from UniProt (or NCBI Protein, where UniProt lacks it) and checks each record. It then writes `data/raw/` and `data/curated/`. `python -m src.verify_data` re-checks the files offline; with `--live` it also checks them against the current UniProt records. The dataset (55 luciferases from 8 families, 12 non-luminous controls), the selection rules and the limitations are documented in [`data/README.md`](data/README.md).

For a full hand-over of the project (decisions, numbers, pitfalls), read [`CONTEXT.md`](CONTEXT.md).

For a plain-language account of methods and results so far, read [`docs/PROJECT_REPORT.md`](docs/PROJECT_REPORT.md).

See [`PLAN.md`](PLAN.md) for scope, syllabus mapping, and execution phases. Background context is in [`docs/bioinformatics_project_handoff.md`](docs/bioinformatics_project_handoff.md).

## Course concepts used

| Week | Topic | Where it appears |
|---|---|---|
| W2, W12 | Motif search, multiple pattern matching | `src/motifs.py`, `results/motif_hits.csv` |
| W3 | Biological databases, GO | UniProt downloads, `manifest.csv`, `results/go_by_family.csv` |
| W5, W14 | Mass spectrometry, peptides | `src/digest.py`, `results/tryptic_summary.csv` |
| W6 | Similarity / homology | `src/align.py`, `results/figures/identity_heatmap.png` |
| W7–W8 | Global/local alignment, DP | `src/align.py`, `results/example_alignments.txt`, `results/dp_example.txt` |
| W9 | DotPlots | `src/dotplot.py`, `results/figures/dotplots.png` |
| W10 | Phylogeny | `src/phylo.py`, `results/figures/trees.png` |
| W11 | Clustering | `src/analysis.py`, `results/figures/pca_features.png` |
| W12 | BLAST | `src/blast.py`, `results/blast_summary.csv` |
| W14 | Proteomics | `src/features.py`, `results/features.csv` |
| W15 | Deep learning (AlphaFold) | `src/structure.py`, `results/alphafold_confidence.csv` |

W1–W2 (intro/motifs via InterPro), W13 (HMMs as Pfam/InterPro profiles), and W15 (deep learning) are mentioned in presentation notes only.

## Stack

Python, Biopython, Pandas, NumPy, Matplotlib, SciPy, scikit-learn, BLAST+

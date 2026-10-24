# Comparative Sequence Analysis of Independently Evolved Luciferases

Course project in bioinformatics: compare luciferase proteins from fireflies, *Renilla*, luminous bacteria, and fungi to ask whether independently evolved light-producing proteins share common sequence-level or physicochemical features.

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
pip install -r requirements.txt
```

See [`PROJECT_PLAN.md`](PROJECT_PLAN.md) for scope, syllabus mapping, and execution phases. Background context is in [`bioinformatics_project_handoff.md`](bioinformatics_project_handoff.md).

## Stack

Python, Biopython, Pandas, NumPy, Matplotlib, SciPy, Seaborn

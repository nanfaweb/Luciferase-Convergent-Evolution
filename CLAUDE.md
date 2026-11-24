# Luciferase convergent evolution

Compares luminous proteins from 10 families (11 independent origins) (plus non-luminous controls) for shared sequence-level or physicochemical features.

Read `PLAN.md` first.

## Stack

Python 3.13: Biopython, pandas, NumPy, SciPy, Matplotlib, scikit-learn. Local BLAST+ 2.17.0 (`~/tools/ncbi-blast-2.17.0+/bin`, on the user PATH).

## Commands

Run from the project root, in this order:

- Data: `python -m src.download`, then `python -m src.verify_data --live`
- Properties: `python -m src.features`
- Similarity: `python -m src.align`
- Dot plots: `python -m src.dotplot`
- BLAST: `python -m src.blast`
- Trees: `python -m src.phylo` (needs `src.align` output)
- Clustering and stats: `python -m src.analysis` (needs `src.features` output)
- Motifs + GO: `python -m src.motifs`; trypsin digest: `python -m src.digest`; AlphaFold: `python -m src.structure`
- Novelty tests: `python -m src.convergence` (pre-registered test, substrate test); `python -m src.pockets` (AlphaFold burial; caches models in `data/structures/`)
- PDF report (after all modules): `python -m src.report_pdf` → `docs/FINAL_REPORT.pdf`
- Notebook: `jupyter nbconvert --to notebook --execute --inplace notebooks/01_analysis.ipynb`

## Project notes

- `data/seeds.csv` is the source of truth for the dataset. Analyse `data/curated/`, never `data/raw/`.
- BLAST+ cannot handle the space in this folder's path, so `src/blast.py` runs it in a system temp folder.
- Keep `docs/PROJECT_REPORT.md` updated after each phase, in simple language.
- Conclusions stay cautious: no claim that a property causes light production.

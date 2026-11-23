# Results

## Protein properties

Rebuild with `python -m src.features`. Values are calculated by Biopython
`ProtParam` from `data/curated/all_sequences.fasta`.

| File | Contents |
|---|---|
| `features.csv` | One row per protein (58), with the columns below |
| `features_by_family.csv` | n, median, min and max of each property, per family and role |

| Column | Meaning |
|---|---|
| `accession`, `family`, `role`, `organism`, `label` | Identity, from the manifest |
| `reviewed` | True for an expert-reviewed Swiss-Prot entry |
| `partial_accepted` | True if the sequence may lack a few end residues, so length and weight may be slightly low |
| `length` | Number of amino acids |
| `molecular_weight_kda` | Molecular weight in kilodaltons |
| `isoelectric_point` | pH at which the protein has no net charge (pI) |
| `gravy` | Average hydrophobicity; negative means more water-loving |
| `aromaticity` | Fraction of residues that are F, W or Y |
| `aa_X_pct` | Percentage of the protein made of amino acid X |

## Reading these values

- **Precursors:** sequences are as deposited, so secreted luciferases
  (copepod, ostracod) still include their signal peptide. This affects
  length, weight and composition slightly.
- **Size:** length and weight mostly reflect protein architecture
  (dinoflagellate luciferases have three repeated domains). GRAVY, pI,
  aromaticity and composition are the fairer cross-family comparisons.
- **Small groups:** `renilla_type` and `ostracod` have 3 luciferases each,
  and some controls are a single protein. Medians here describe the data;
  they are not statistical tests.

## Analysis outputs

| File | Made by | Contents |
|---|---|---|
| `pairwise.csv` | `src.align` | Every protein pair: global identity %, local score, z-score vs shuffles |
| `identity_matrix.csv` | `src.align` | 58 x 58 global identity matrix |
| `similarity_summary.csv` | `src.align` | Identity and z-score by pair type |
| `example_alignments.txt` | `src.align` | NW and SW alignments for a same-family pair and a cross-family pair |
| `blast_all_vs_all.tsv` | `src.blast` | All BLAST hits at E < 0.001 |
| `blast_summary.csv` | `src.blast` | Same-family and cross-family hit counts per protein |
| `trees/<family>.nwk` | `src.phylo` | Within-family trees (Newick format) |
| `stats.csv` | `src.analysis` | Clustering ARI, Kruskal-Wallis and sign tests |
| `control_comparison.csv` | `src.analysis` | Luciferase vs control medians per family |
| `figures/` | all | Heatmap, dot plots, trees, PCA and boxplots |

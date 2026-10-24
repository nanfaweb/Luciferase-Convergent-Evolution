# Final Plan: Comparative Sequence Analysis of Independently Evolved Luciferases

## 1. Project identity

| Item | Decision |
|---|---|
| Full title | Molecular Signatures of Convergent Evolution: Comparative Sequence Analysis of Independently Evolved Luciferases |
| Question | Do independently evolved luciferases share common sequence-level or physicochemical characteristics despite different evolutionary origins? |
| Type | Undergraduate CS elective project/presentation — not a publication |
| Core method | Protein-sequence comparative analysis in Python (not DNA motifs, not ML, not wet lab) |
| Proposal | Already exists — do not rewrite; use as scope contract |
| Source of truth | [bioinformatics_project_handoff.md](bioinformatics_project_handoff.md) |

**Plain-language pitch:** Compare luciferase proteins from fireflies, marine animals (*Renilla*), luminous bacteria, and fungi using public sequences, alignments, BLAST, DotPlots, within-family trees, and measurable protein properties to see what is shared vs family-specific.

## 2. Biological groups and data

**Luciferase groups (target ~5–10 sequences each after expansion):**
1. Firefly luciferase
2. *Renilla* / sea pansy luciferase
3. Bacterial LuxA (alpha subunit)
4. Fungal Luz luciferase

**Starter accessions (must re-validate on UniProt before locking the dataset):**

- Firefly: P08659, P13129, Q5UFR2
- Renilla: P27652, Q9BLZ3
- Bacterial LuxA: P19907, P09140
- Fungal Luz: A0A3G9JYH7, A0A2H3E985
- Controls (non-luminous): Q9VCC6 (acyl-CoA synthetase), P22643 (haloalkane dehalogenase DhlA)

**Download pattern:** `https://rest.uniprot.org/uniprotkb/{ACCESSION}.fasta`

**Scientific constraints (always):**
- Many luciferase systems evolved independently; do not claim every pair is wholly unrelated.
- Firefly luciferases relate to acyl-CoA synthetases; *Renilla* relates to haloalkane-dehalogenase-like recruitment (parallel evolution).
- No single MSA/phylogenetic tree that forces all four families together.
- Low pairwise identity does not prove absence of homology; controls are suggestive, not causal proof of light-related properties.

## 3. Repository structure

```
Bio Proj/
  data/
    raw/                 # per-accession FASTA downloads
    curated/
      manifest.csv       # family, species, accession, URLs, length, GO/InterPro, reason
      luciferases.fasta  # combined curated sequences
      controls.fasta
  notebooks/
    01_analysis.ipynb    # instructor-facing walkthrough
  src/
    download.py
    features.py
    align.py
    dotplot.py
    blast_utils.py
    phylo_within.py
    cluster_features.py
    plots.py
  results/
    features.csv
    similarity_matrix.csv
    figures/             # heatmap, DotPlots, trees, cluster plot
    blast_summary.md
  requirements.txt
  README.md              # setup + "Course concepts used"
```

**Stack:** Python 3, Biopython, Pandas, NumPy, Matplotlib, SciPy, Seaborn. BLAST via Biopython/NCBI web export if BLAST+ is unavailable. No GPU, deep learning, or web app.

## 4. Syllabus coverage (instructor map)

**Implement in the pipeline**

- **W3 Databases / GO** — UniProt/NCBI downloads; manifest; InterPro/GO notes for representatives
- **W6 Similarity / homology** — pairwise identity matrix; within- vs across-family interpretation
- **W7–W8 Global/local alignment + DP** — Needleman–Wunsch and Smith–Waterman on selected pairs; short DP/gap-penalty note in notebook
- **W9 DotPlots** — within-family pair(s) vs one cross-family pair
- **W10 Phylogeny** — within-family MSA + distance tree only
- **W12 BLAST** — BLASTP query against curated set / public DB; hit summary
- **W11 (light) Clustering** — cluster proteins by physicochemical feature vectors
- **W14 Proteomics** — ProtParam: length, MW, AA composition, GRAVY, aromaticity, pI

**Mention in presentation only**

- **W1–W2** Intro / motifs → protein sequences + InterPro domains (not DNA clock motifs)
- **W13** HMMs → Pfam/InterPro profiles underlie domain annotations; we do not train HMMs
- **W15** Deep learning → out of scope; optional “future work” line if asked

**Do not implement**

- W4 genome assembly; W5 antibiotics/MS; W9 sorting-by-reversals; W11 transcriptomics/expression; W13 HIV phenotyping / HMM training; W15 DL/image models; cross-family forced phylogeny

## 5. Pipeline (execution order)

```mermaid
flowchart TD
  phase1[Phase1_Scaffold_and_Data]
  phase2[Phase2_Core_Analyses]
  phase3[Phase3_Expand_and_Polish]
  phase1 --> phase2 --> phase3

  subgraph p1 [Phase1]
    scaffold[Scaffold_repo]
    validate[Validate_accessions]
    manifest[Write_manifest]
    download[Download_FASTA]
  end

  subgraph p2 [Phase2]
    features[ProtParam_features]
    align[Global_Local_alignment]
    heatmap[Similarity_heatmap]
    dots[DotPlots]
    blast[BLASTP]
    phylo[Within_family_trees]
    cluster[Feature_clustering]
  end

  subgraph p3 [Phase3]
    expand[Expand_to_5to10_per_group]
    notebook[Polish_notebook_narrative]
    slides[Presentation_syllabus_slide]
    conclude[Cautious_conclusions]
  end
```

### Phase 1 — Scaffold and starter data
1. Create folders, `requirements.txt`, README with syllabus map.
2. Live-check each starter + control UniProt record (name, organism, completeness).
3. Write `manifest.csv`; download FASTA to `data/raw/`; build curated FASTA files.

### Phase 2 — Core analyses on starter set
1. **Features:** Biopython `ProtParam` → `results/features.csv` + family comparison plots.
2. **Alignment:** pairwise global/local for selected pairs; full pairwise identity → heatmap.
3. **DotPlots:** 1–2 within-family + 1 cross-family figure.
4. **BLAST:** one luciferase query; document top hits and family membership.
5. **Phylogeny:** within-family tree(s) only where ≥3–4 sequences exist (expand first if needed).
6. **Clustering:** hierarchical or k-means on feature vectors; color by family.

### Phase 3 — Expand dataset and course deliverables
1. Grow each luciferase group toward ~5–10 sequences; re-run analyses.
2. Finish `01_analysis.ipynb` in course language (alignment, BLAST, DotPlot, phylogeny, clustering).
3. Package figures + short limitations/conclusions; one presentation slide mapping weeks → outputs.

## 6. Expected deliverables checklist

- Curated FASTA + documented `manifest.csv` (accessions, sources)
- Reproducible Python scripts + analysis notebook
- Pairwise similarity matrix and heatmap
- DotPlot figures; BLAST summary
- Within-family phylogenetic tree(s)
- Physicochemical feature table, plots, and feature clustering figure
- Control comparisons labeled as controls
- Short interpretation with limitations; presentation materials including syllabus coverage slide

## 7. Work checklist

1. Scaffold folders, `requirements.txt`, README with syllabus map
2. Validate UniProt accessions; write `manifest.csv` with GO/InterPro notes
3. Download and curate starter FASTA dataset
4. ProtParam features table + feature-space clustering plots
5. Global/local pairwise alignment, identity heatmap, DotPlots
6. BLASTP summary + within-family MSA/distance trees
7. Expand to ~5–10 sequences per luciferase family + controls
8. Notebook narrative, figures, cautious conclusions, presentation syllabus slide

## 8. Interpretation rules (for results write-up)

- Expect high within-family similarity and low cross-family sequence similarity.
- Physicochemical “shared patterns,” if any, are hypotheses — do not claim they cause bioluminescence.
- Phrase carefully: “Many luciferase systems evolved independently…”; analysis asks whether families remain sequence-distinct while checking for any common property patterns.

## 9. Immediate first actions

1. Scaffold the repository layout and dependencies.
2. Validate the 9 starter + 2 control accessions and write the first manifest.
3. Download those FASTA files and smoke-test features + one alignment + one DotPlot before dataset expansion.

## 10. Out of scope

DNA motif / molecular-clock search as a core aim; genome assembly; transcriptomics; neural models; wet lab; paper publication; one global tree of all luciferase families.

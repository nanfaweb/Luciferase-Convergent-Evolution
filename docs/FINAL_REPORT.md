# Molecular Signatures of Convergent Evolution: Comparative Sequence Analysis of Independently Evolved Luciferases

**Final report.** Bioinformatics course project, FAST-NU Lahore, 2026.
Code, data and figures are in this repository. Re-run everything with
`notebooks/01_analysis.ipynb`. A PDF version with all figures is
`docs/FINAL_REPORT.pdf`.

---

## 1. Problem formulation

**Background.** Bioluminescence, the production of light by living organisms,
has evolved independently many times. In each case a light-emitting protein,
usually a **luciferase**, oxidises a substrate (a luciferin) to make light.
Different lineages recruited these proteins from unrelated protein families:
- acyl-CoA synthetases (ACS) in beetles;
- haloalkane dehalogenases (HLD) in the sea pansy *Renilla*;
- flavin-dependent monooxygenases in bacteria;
- calcium-binding EF-hand proteins in jellyfish photoproteins;
- lipocalin-like proteins in a deep-sea shrimp.

These proteins perform the same *function*. Whether they also converged at
the *molecular* level is an open, testable question.

**Research question.** Do independently evolved luminous proteins share
protein-sequence or physicochemical characteristics despite their different
evolutionary origins?

**Hypotheses.**
- **H1 (sequence).** Luminous proteins from independent origins share no
  detectable sequence similarity, motifs or peptides.
- **H2 (properties).** If light production imposes shared constraints, each
  luminous protein should differ from its closest non-luminous relatives
  (controls) **in the same direction** across independent families.
- **H0 for H2.** Properties are set by protein-family ancestry, and the
  luminous-vs-control differences point in random directions.

**Scope.** This is a protein-level computational comparison of public
sequences and predicted structures. There is no wet-lab work, and no claim
that any property *causes* light emission.

---

## 2. Existing solutions and related work

**Luciferase evolution**
- Delroisse et al. [1] reviewed luciferases and photoproteins across phyla and
  described many independent evolutionary origins.
- He et al. [2] showed that fireflies and click beetles each recruited
  luciferase from the ACS family. Zhang et al. [21] found that luminous
  beetle luciferases share a distinct residue pattern at the luciferin-binding
  site.
- Lau et al. [4] showed that a brittle star recruited an HLD as a luciferase
  independently of *Renilla*.
- Other systems are each characterised in their own study:
  - fungal Luz, Kotlobay et al. [3];
  - copepods, Takenaka et al. [5];
  - dinoflagellates, Liu et al. [6];
  - Ca²⁺-regulated photoproteins such as aequorin [23, 24];
  - the Oplophorus shrimp luciferase [22].

**Standard methods**
- Pairwise alignment by dynamic programming, global [7] and local [8], scored
  with BLOSUM62 [9].
- BLAST+ [10] for homology search.
- Neighbour joining [11] for distance trees.
- ProtParam [13] and Kyte–Doolittle hydropathy [14] for physicochemical
  properties.
- PROSITE [18] for motifs.
- Pfam/InterPro profile HMMs for protein domains.
- AlphaFold [15, 16] for predicted structures.

**What is missing.** Earlier work studies **one family at a time**, or compares
families **qualitatively**: different folds, different luciferins.
- Residue-level convergence has been described *within* beetles [2, 21].
- In a literature search on 2026-11-24 we did **not** find a quantitative,
  cross-family test of whether luminous proteins differ from their own
  non-luminous relatives in a consistent physicochemical direction.
- We also did not find a test of whether unrelated proteins that use the
  *same* luciferin (coelenterazine) converge in composition.

That search was limited, so we say "did not find", not "first".

---

## 3. Proposed methodology (framework)

```
UniProt / NCBI ──► curation (seeds.csv rules) ──► 81 proteins: 62 luminous in 10 families + 19 controls
                                                   │
        ┌──────────────────┬───────────────────────┼──────────────────────┬────────────────────┐
        ▼                  ▼                       ▼                      ▼                    ▼
  Sequence level     Homology search         Motifs & peptides     Properties (ProtParam)   Structure
  NW / SW + shuffle  BLAST+ all-vs-all       PROSITE multi-pattern GRAVY, pI, aromaticity,  AlphaFold DB:
  z-scores, dot plots                        scan; trypsin digest  composition              confidence, burial
  within-family NJ trees                                                 │
                                                                         ▼
                 clustering (PCA, Ward + ARI) · Kruskal–Wallis · luminous-vs-own-control tests
                 NOVEL: pre-registered replication · within-family permutation test ·
                        substrate-convergence test · surface-vs-core localisation
```

**Key design choices**
1. **Matched controls.** Each of 7 families has 1–4 non-luminous homologs from
   the same protein family or superfamily. Comparing a luminous protein with
   its own controls cancels out family-level traits.
2. **The family is the unit of evidence.** Sign tests count families, and
   permutation tests shuffle labels only *within* families. This way, shared
   ancestry cannot create a signal.
3. **Pre-registration.** After 5 families suggested a pattern, the prediction
   and the pass/fail rules were committed to git (commit `d2978e2`) before
   any data for the new families was downloaded. The controls were chosen by
   a fixed rule (reviewed entries, same superfamily, one per phylum or order)
   without looking at their properties.
4. **No mixed-family tree.** Unrelated proteins cannot be meaningfully aligned
   together.
5. **Composition, not size.** Cross-family comparisons use composition-type
   features.

---

## 4. Implementation and experiments

**Environment.**
- Python 3.13 with Biopython 1.86 [12], pandas, NumPy, SciPy, Matplotlib and
  scikit-learn.
- Local BLAST+ 2.17.0.
- Data from the UniProt [17], NCBI, PROSITE and AlphaFold DB REST APIs.

### 4.1 Dataset (`src/download.py`, `src/verify_data.py`)

**Selection rules.**
- There must be functional evidence that the protein emits light.
- The sequence must be complete. 6 near-complete clones and 3 mature-chain
  controls are kept and flagged.
- Within a family, no two sequences may be ≥97% identical.
- Expert-reviewed entries are preferred.
- Every record is versioned, hashed and re-verified against live UniProt.

| Family | Luminous | Controls | Origin |
|---|---|---|---|
| Firefly | 8 | 2 (*Drosophila* CG6178, *Tenebrio* LL) | ACS family |
| Click beetle | 4 | 1 (non-luminous click beetle AbLL) | ACS family, independent of fireflies [2] |
| *Renilla*-type | 3 (2 *Renilla*, 1 brittle star) | 3 (HLDs LinB, DhaA, DafA) | HLD family, two independent recruitments [4] |
| Bacterial LuxA | 10 | 3 (SsuD) | Luciferase-like monooxygenase family |
| Fungal Luz | 9 | 3 (Luz-domain proteins, non-luminous fungi) | Luz family [3] |
| Photoprotein (new) | 6 (aequorin ×2, mitrocomin, obelin ×2, clytin) | 4 (3 sarcoplasmic Ca-binding proteins, calmodulin) | EF-hand family [23, 24] |
| Oplophorus shrimp (new) | 1 | 3 (fatty-acid-binding proteins) | Calycin / lipocalin superfamily [22] |
| Copepod | 10 | – | Unique family [5] |
| Dinoflagellate | 8 | – | Unique three-domain enzyme [6] |
| Ostracod | 3 | – | Cypridinid luciferase |
| **Total** | **62** | **19** | **11 independent origins** |

Photoproteins are Ca²⁺-triggered and not classic turnover luciferases, so
every result is also reported without them.

### 4.2 Experiments

| # | Experiment | Module | Settings |
|---|---|---|---|
| E1 | Physicochemical properties | `features.py` | ProtParam: length, MW, pI, GRAVY, aromaticity, amino-acid % |
| E2 | Pairwise alignment, all 3,240 pairs | `align.py` | NW identity; SW score with z-score vs 20 shuffles; BLOSUM62, gap −10/−0.5 |
| E3 | Dynamic-programming worked example | `align.py` | Hand-written NW table, checked against Biopython |
| E4 | Dot plots | `dotplot.py` | 10-residue window, BLOSUM62 ≥ 23 |
| E5 | BLAST all-vs-all | `blast.py` | blastp, E < 10⁻³ |
| E6 | Within-family trees | `phylo.py` | Poisson distance, neighbour joining, rooted on a control |
| E7 | Motif search (multiple pattern matching) | `motifs.py` | All PROSITE patterns UniProt lists for the dataset, scanned as regexes |
| E8 | In-silico trypsin digest | `digest.py` | Cut after K/R unless the next residue is P; 6–30 aa; monoisotopic mass |
| E9 | Clustering and statistics | `analysis.py` | PCA; Ward vs family (ARI); Kruskal–Wallis; sign test |
| E10 | AlphaFold confidence | `structure.py` | Mean pLDDT from AlphaFold DB |
| **E11** | **Pre-registered convergence test** | `convergence.py` | One-sided sign test; within-family label permutation (10,000, seed 1); Holm; robustness; calibration on 500 random features |
| **E12** | **Substrate convergence** (exploratory) | `convergence.py` | Family centroids in standardized composition space; exact test over all 210 label sets |
| **E13** | **Surface vs core** (exploratory) | `pockets.py` | AlphaFold models; Cα neighbours within 10 Å; pLDDT ≥ 70; per-protein median split; two-sided permutation |

---

## 5. Results

### 5.1 Sequence similarity (E2, E4, E5)

| Pair type | Pairs | Median identity | Median SW z-score |
|---|---|---|---|
| Same family | 209 | 68% | 130 |
| Luminous protein vs its own control | 113 | 30% | 15 |
| Firefly vs click beetle (same ACS superfamily) | 32 | 49% | 137 |
| All other cross-family luminous pairs | 1,650 | 20% | 0.0 (chance) |

- **BLAST:** 684 significant hits. The only cross-family hits link fireflies
  and click beetles.
- **Weak alignment hits:** 45 unrelated pairs reached z > 3 (maximum 6.7).
  This is expected by chance across more than 1,600 tests, and BLAST
  supported none of them.
- **Dot plots:** a continuous diagonal within a family (479 dots), a broken
  one for luciferase vs control (334), and noise for unrelated proteins (6).

**E3.** The hand-written NW table scores −8, the same as Biopython.

**E6.** The trees reproduce known biology. The brittle-star luciferase groups
with its own non-luminous dehalogenase, not with *Renilla*.

**E7.** The dataset has 4 PROSITE patterns. None crosses independent
families:
- AMP-binding: firefly and click beetle, including controls;
- bacterial luciferase: LuxA only, not SsuD;
- EF-hand: photoproteins and their EF-hand controls;
- FABP: only the FABP controls. The Oplophorus luciferase is too divergent
  to match it.

**E8.** Of 1,976 tryptic peptides, only 2 occur in more than one family, and
both are shared by fireflies and click beetles.

**E1/E9.**
- Composition recovers the families: ARI = 0.82, where 1 is perfect.
- Every main property differs strongly between families: Kruskal–Wallis
  p < 10⁻⁶.

**E10.** 76 of 81 proteins have AlphaFold models.
- Median pLDDT is about 90 for most families.
- Oplophorus: 85. Dinoflagellates: 86. Copepods: 56.

### 5.2 Novel analysis 1: pre-registered replication test (E11)

**Prediction, written before the new data:** each new luminous protein has
lower GRAVY and higher aromaticity than its controls.

| New family | GRAVY (predicted lower) | Aromaticity (predicted higher) |
|---|---|---|
| Photoprotein | lower ✔ (−0.32) | lower ✘ (−0.007) |
| Oplophorus shrimp | **higher ✘ (+0.76)** | higher ✔ (+0.032) |

**The pre-registered replication failed.** Each new family contradicted one
of the two predictions.

Main test across all 7 families with controls:

| Analysis | Property | Families in predicted direction | Sign test p | Permutation p (Holm) |
|---|---|---|---|---|
| All 7 families | GRAVY lower | 6/7 | 0.0625 | 0.55 (0.55) |
| All 7 families | **Aromaticity higher** | **6/7** | 0.0625 | **0.010 (0.021)** |
| Firefly + click beetle as one origin | Aromaticity | 5/6 | 0.11 | 0.014 (0.028) |
| Photoproteins excluded | Aromaticity | 6/6 | 0.016 | 0.009 (0.019) |
| Leave-one-family-out (7 runs) | Aromaticity | 5/6 or 6/6 | – | 0.009–0.072 |
| Original 5 families only (the hint) | GRAVY / aromaticity | 5/5 / 5/5 | 0.031 / 0.031 | 0.0001 / 0.0003 |

- **The test behaves correctly.** On 500 random features, 4% gave p < 0.05,
  close to the 5% expected.
- **GRAVY is not supported.** Six families point the predicted way, but the
  Oplophorus luciferase is much *more* hydrophobic than its controls
  (+3.7 SD), and the test result is p = 0.55.
- **Aromaticity is supported, with caveats.**
  - It is significant after Holm correction, and robust to merging the two
    beetle families and to removing photoproteins.
  - It weakens when LuxA is left out (p = 0.072).
  - Of the two new families, only Oplophorus agreed.
- **Cysteine %** (not predicted, two-sided) also shifts: p = 0.013.

### 5.3 Novel analysis 2: same substrate, same composition? (E12, exploratory)

Four unrelated folds use the same luciferin, coelenterazine: *Renilla*-type,
copepod, photoprotein and Oplophorus. If the substrate shaped composition,
these families should sit closer together than random sets of 4 families.

They sit **further apart**: mean centroid distance 8.2, against a null
median of 6.3; exact p = 0.995 over 210 label sets.

**Sharing a luciferin does not make unrelated proteins compositionally
similar.**

### 5.4 Novel analysis 3: surface or core? (E13, exploratory)

Each residue of each AlphaFold model was classed as buried or exposed. The
aromatic excess of luminous proteins over their controls is:
- **on exposed residues:** higher in 6 of 7 families, including both new
  families (photoprotein +4.4 points, Oplophorus +6.4 points). Permutation
  p = 0.0001.
- **not on buried residues:** p = 0.70.
- Hydropathy shows no buried or exposed difference (p = 0.98 and 0.45).

So the photoproteins, which failed the whole-protein aromaticity prediction,
*do* carry more aromatic residues on their surface. They lose them in the
core.

Caveats:
- "Exposed" includes the surface-facing linings of binding clefts, so this
  proxy cannot separate pockets from the general surface.
- This analysis was designed after the E11 results, so it needs independent
  confirmation.

Figures in `results/figures/`:
- `similarity_by_category.png`, `identity_heatmap.png`, `dotplots.png`;
- `trees.png`, `pca_features.png`, `feature_boxplots.png`;
- `control_differences.png`, `alphafold_plddt.png`;
- `h2_replication.png`, `substrate_distance.png`, `aromatic_burial.png`.

---

## 6. Conclusion

1. **H1 is supported: no shared sequence.** Across 11 independent origins,
   luminous proteins share no detectable sequence, motifs or peptides. The
   only exception is firefly/click beetle, which share an ACS ancestor.
2. **Properties follow ancestry.** Composition sorts proteins by family
   (ARI 0.82).
3. **New: a pre-registered test rejected the "less hydrophobic" hint.** It
   looked consistent in 5 families but did not generalise. This is a
   negative result, and reporting it is part of the contribution.
4. **New: a family-controlled aromaticity association.** Across 7 families,
   luminous proteins are more aromatic than their own non-luminous relatives
   (p = 0.010, Holm 0.021). It replicated in one new family (Oplophorus) but
   not the other (photoproteins).
5. **New (exploratory): the excess is on the surface.** In 3D models, the
   aromatic excess sits on exposed residues (p = 0.0001), including in both
   new families, and not in the core.
6. **New: sharing a substrate does not cause compositional convergence.**
   Coelenterazine users from four unrelated folds are no more similar than
   random families.

**Answer.** Independently evolved luminous proteins converged in function,
not in sequence. Their composition is set mainly by ancestry. There is one
cautious exception: a modest, family-controlled enrichment of aromatic
residues, located mostly on the protein surface. This is an
**association**, not a cause. Tyr and Trp are common in luciferins and in
protein–ligand interfaces, which suggests a hypothesis about substrate or
partner binding that would need experiments to test.

**Contribution beyond the literature**
- The first test (that we found) of the cross-family "luminous vs own
  non-luminous relative" property question with family-stratified
  statistics.
- A pre-registered replication, honestly reported as failed for GRAVY.
- A substrate-convergence test, answered negatively.
- A structural localisation of the aromatic signal.

**Limitations**
- 7 families with controls is still a small number. One family (Oplophorus)
  has a single luminous sequence, and its FABP controls are only distantly
  related.
- Photoproteins are not classic luciferases.
- Secreted proteins are analysed with their signal peptides.
- The burial measure is a proxy, and E12 and E13 are exploratory.

**Future work**
- Confirm the surface-aromatic result on further families with controls,
  such as ctenophore photoproteins, *Odontosyllis* or the tunicate
  luciferase.
- Use real pocket detection and compare the residues that line the
  luciferin-binding pocket.
- Analyse mature proteins, with signal peptides removed.

### Course coverage

| Week | Topic | Used in this project |
|---|---|---|
| 1 | Intro, origin of replication | Background only (genome-level; not applicable) |
| 2 | Motif search | PROSITE motif scan (E7) |
| 3 | Biological databases, GO | UniProt/NCBI curation, InterPro/Pfam, GO summary |
| 4 | Genome assembly | Not applicable (no genome data) |
| 5 | Mass spectrometry | Trypsin digest and peptide masses (E8) |
| 6 | Similarity and homology | Identity, z-scores and BLAST (E2, E5) |
| 7 | Global and local alignment | Needleman–Wunsch and Smith–Waterman (E2) |
| 8 | Dynamic programming | Hand-written DP table and traceback (E3) |
| 9 | Dot plots | Dot plots (E4) |
| 10 | Phylogenetics | Neighbour-joining trees (E6) |
| 11 | Clustering | PCA, Ward, ARI (E9); composition distances (E12) |
| 12 | Multiple pattern matching, BLAST | Multi-pattern PROSITE scan (E7), BLAST+ (E5) |
| 13 | Hidden Markov Models | Pfam/InterPro profile-HMM domains for curation and for choosing controls |
| 14 | Proteomics, peptide sequencing | ProtParam (E1), tryptic peptides (E8) |
| 15 | Deep learning for biomolecules | AlphaFold confidence (E10) and structure-based burial analysis (E13) |

---

## References

1. Delroisse J, Duchatelet L, Flammang P, Mallefet J. Leaving the dark side? Insights into the evolution of luciferases. *Frontiers in Marine Science* 8:673620 (2021). https://doi.org/10.3389/fmars.2021.673620
2. He J, et al. Multiple origins of bioluminescence in beetles and evolution of luciferase function. *Molecular Biology and Evolution* 41(1):msad287 (2024). https://doi.org/10.1093/molbev/msad287
3. Kotlobay AA, et al. Genetically encodable bioluminescent system from fungi. *PNAS* 115(50):12728–12732 (2018). https://doi.org/10.1073/pnas.1803615115
4. Lau ES, et al. Functional characterization of luciferase in a brittle star indicates parallel evolution influenced by genomic availability of haloalkane dehalogenase. *Molecular Biology and Evolution* 42(5):msaf081 (2025). https://doi.org/10.1093/molbev/msaf081
5. Takenaka Y, et al. Evolution of bioluminescence in marine planktonic copepods. *Molecular Biology and Evolution* 29(6):1669–1681 (2012). https://doi.org/10.1093/molbev/mss009
6. Liu L, Wilson T, Hastings JW. Molecular evolution of dinoflagellate luciferases, enzymes with three catalytic domains in a single polypeptide. *PNAS* 101(47):16555–16560 (2004). https://doi.org/10.1073/pnas.0407597101
7. Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology* 48(3):443–453 (1970). https://doi.org/10.1016/0022-2836(70)90057-4
8. Smith TF, Waterman MS. Identification of common molecular subsequences. *Journal of Molecular Biology* 147(1):195–197 (1981). https://doi.org/10.1016/0022-2836(81)90087-5
9. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. *PNAS* 89(22):10915–10919 (1992). https://doi.org/10.1073/pnas.89.22.10915
10. Camacho C, et al. BLAST+: architecture and applications. *BMC Bioinformatics* 10:421 (2009). https://doi.org/10.1186/1471-2105-10-421
11. Saitou N, Nei M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4(4):406–425 (1987). https://doi.org/10.1093/oxfordjournals.molbev.a040454
12. Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423 (2009). https://doi.org/10.1093/bioinformatics/btp163
13. Gasteiger E, et al. Protein identification and analysis tools on the ExPASy server. In: Walker JM (ed.), *The Proteomics Protocols Handbook*, Humana Press, 571–607 (2005). https://doi.org/10.1385/1-59259-890-0:571
14. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157(1):105–132 (1982). https://doi.org/10.1016/0022-2836(82)90515-0
15. Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2
16. Varadi M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Research* 50(D1):D439–D444 (2022). https://doi.org/10.1093/nar/gkab1061
17. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51(D1):D523–D531 (2023). https://doi.org/10.1093/nar/gkac1052
18. Sigrist CJA, et al. New and continuing developments at PROSITE. *Nucleic Acids Research* 41(D1):D344–D347 (2013). https://doi.org/10.1093/nar/gks1067
19. Kruskal WH, Wallis WA. Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association* 47(260):583–621 (1952). https://doi.org/10.1080/01621459.1952.10483441
20. Hubert L, Arabie P. Comparing partitions. *Journal of Classification* 2(1):193–218 (1985). https://doi.org/10.1007/BF01908075

21. Zhang R, He J, Dong Z, et al. Genomic and experimental data provide new insights into luciferin biosynthesis and bioluminescence evolution in fireflies. *Scientific Reports* 10:15882 (2020). https://doi.org/10.1038/s41598-020-72900-z
22. Inouye S, et al. Secretional luciferase of the luminous shrimp *Oplophorus gracilirostris*: cDNA cloning of a novel imidazopyrazinone luciferase. *FEBS Letters* 481(1):19–25 (2000). https://doi.org/10.1016/S0014-5793(00)01963-3
23. Inouye S, et al. Cloning and sequence analysis of cDNA for the luminescent protein aequorin. *PNAS* 82(10):3154–3158 (1985). https://doi.org/10.1073/pnas.82.10.3154
24. Shimomura O, Johnson FH, Saiga Y. Extraction, purification and properties of aequorin, a bioluminescent protein from the luminous hydromedusan, *Aequorea*. *Journal of Cellular and Comparative Physiology* 59(3):223–239 (1962). https://doi.org/10.1002/jcp.1030590302

*All 24 references were checked (1-20 on 2026-11-23, 21-24 on 2026-11-24) against their Crossref records
(authors, title, journal, volume, pages, year, DOI).*

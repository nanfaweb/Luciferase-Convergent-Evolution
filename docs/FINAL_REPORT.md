# Molecular Signatures of Convergent Evolution: Comparative Sequence Analysis of Independently Evolved Luciferases

**Final report.** Bioinformatics course project, FAST-NU Lahore, 2026.
Code, data and figures: this repository. Re-run everything with `notebooks/01_analysis.ipynb`.

---

## 1. Problem formulation

**Background.** Bioluminescence, the production of light by living organisms,
has evolved independently many times: in beetles, marine animals, bacteria,
fungi, crustaceans and dinoflagellates. In each case an enzyme called a
**luciferase** drives a light-emitting oxidation reaction. Different lineages
recruited luciferases from unrelated protein families: acyl-CoA synthetases
in beetles, haloalkane dehalogenases in the sea pansy *Renilla*, and
flavin-dependent monooxygenases in bacteria. These enzymes perform the same
*function*. Whether they also converged at the *molecular* level is an open,
testable question.

**Research question.** Do independently evolved luciferases share common
protein-sequence or physicochemical characteristics despite their different
evolutionary origins?

**Hypotheses.**
- **H1 (sequence).** Luciferases from independent origins share no detectable
  sequence similarity, motifs or peptides.
- **H2 (properties).** If light production imposes shared constraints,
  luciferases should differ from their closest non-luminous relatives
  (controls) **in the same direction** across independent families.
- **H0 for H2.** Properties are set by protein-family ancestry. Luciferases
  group by family, and the luciferase-vs-control differences point in random
  directions.

**Scope.** This is a protein-level computational comparison of public
sequences. It involves no DNA motif search, no wet-lab work, and no claim
that any property *causes* light emission.

---

## 2. Existing solutions and related work

- **Luciferase evolution.** Delroisse et al. [1] compared luciferases and
  photoproteins across many phyla and described several independent
  evolutionary groups. He et al. [2] showed that beetle bioluminescence has
  multiple origins: fireflies and click beetles each recruited luciferase
  from the acyl-CoA synthetase (ACS) family. Lau et al. [4] showed that a
  brittle star recruited a haloalkane dehalogenase as a luciferase
  independently of *Renilla*, a case of parallel evolution. Kotlobay et al.
  [3] identified the fungal luciferase Luz as a distinct protein family.
  Takenaka et al. [5] and Liu et al. [6] characterised copepod and
  dinoflagellate luciferases, two further unrelated systems.
- **Standard methods.**
  - **Pairwise alignment:** dynamic programming, global [7] and local [8],
    with BLOSUM62 scoring [9].
  - **Homology search:** BLAST+ [10].
  - **Distance trees:** neighbour joining [11].
  - **Physicochemical profiling:** ExPASy ProtParam [13] and Kyte–Doolittle
    hydropathy [14].
  - **Motif databases:** PROSITE [18].
  - **Domain annotation:** profile HMMs, through Pfam and InterPro.
  - **Deep-learning structure prediction:** AlphaFold [15, 16].
- **The gap this project addresses.** Previous studies mostly built family
  trees or characterised single enzymes. We combine every standard method on
  **one curated, multi-origin dataset with matched non-luminous controls**.
  That makes it possible to separate "inherited from the protein family" from
  "shared because the protein makes light".

---

## 3. Proposed methodology (framework)

```
UniProt / NCBI ──► curation (seeds.csv rules) ──► 67 proteins: 55 luciferases in 8 families + 12 controls
                                                   │
        ┌──────────────────┬───────────────────────┼──────────────────────┬────────────────────┐
        ▼                  ▼                       ▼                      ▼                    ▼
  Sequence level     Homology search         Motifs & peptides     Properties (ProtParam)   Structure
  NW / SW + shuffle  BLAST+ all-vs-all       PROSITE multi-pattern GRAVY, pI, aromaticity,  AlphaFold DB
  z-scores, dot plots                        scan; trypsin digest  composition              confidence
  within-family NJ trees                                                 │
                                                                         ▼
                                    clustering (PCA, Ward + ARI) · Kruskal–Wallis ·
                                    luciferase-vs-own-control sign test (unit = family)
```

**Key design choices.**
1. **Matched controls.** Each of 5 families has 1–3 non-luminous homologs
   from the same protein family. Comparing a luciferase with its own control
   cancels out family-level traits.
2. **The family is the unit of evidence.** Sequences within a family share
   ancestry, so the statistics compare families, not individual proteins.
3. **No mixed-family tree.** Unrelated proteins cannot be meaningfully
   aligned together, so trees are built within families only.
4. **Composition, not size.** Length and weight reflect domain architecture,
   so cross-family comparisons use composition-type features.

---

## 4. Implementation and experiments

**Environment.**
- Python 3.13 with Biopython 1.86 [12], pandas, NumPy, SciPy, Matplotlib and
  scikit-learn.
- Local BLAST+ 2.17.0.
- Data from the UniProt [17], NCBI Protein, PROSITE and AlphaFold DB REST APIs.

### 4.1 Dataset (`src/download.py`, `src/verify_data.py`)

**Selection rules.**
- There must be functional evidence that the protein is a luciferase.
- The sequence must be complete; 6 near-complete clones are kept and flagged.
- Within a family, no two sequences may be ≥97% identical (one per species).
- Expert-reviewed entries are preferred.
- Every record is versioned and hashed, and re-verified against live UniProt.

| Family | Luciferases | Controls | Origin of the luciferase |
|---|---|---|---|
| Firefly | 8 | 2 (*Drosophila* CG6178, *Tenebrio* LL) | ACS family |
| Click beetle | 4 | 1 (non-luminous click beetle AbLL) | ACS family, independent of fireflies [2] |
| *Renilla*-type | 3 (2 *Renilla*, 1 brittle star) | 3 (dehalogenases LinB, DhaA, DafA) | HLD family, two independent recruitments [4] |
| Bacterial LuxA | 10 | 3 (SsuD from 3 species) | Luciferase-like monooxygenase family |
| Fungal Luz | 9 | 3 (Luz-domain proteins from non-luminous fungi) | Luz family [3] |
| Copepod | 10 | – | Unique family [5] |
| Dinoflagellate | 8 | – | Unique three-domain enzyme [6] |
| Ostracod | 3 | – | Cypridinid luciferase |
| **Total** | **55** | **12** | **9 independent origins** |

### 4.2 Experiments

| # | Experiment | Module | Settings |
|---|---|---|---|
| E1 | Physicochemical properties | `features.py` | ProtParam: length, MW, pI, GRAVY, aromaticity, amino-acid % |
| E2 | Pairwise alignment, all 2,211 pairs | `align.py` | Needleman–Wunsch identity; Smith–Waterman score with z-score vs 20 shuffles; BLOSUM62, gap −10/−0.5 |
| E3 | Dynamic-programming worked example | `align.py` | Hand-written NW table for HEAGAWGHEE vs PAWHEAE, checked against Biopython |
| E4 | Dot plots | `dotplot.py` | 10-residue window, BLOSUM62 ≥ 23 |
| E5 | BLAST all-vs-all | `blast.py` | blastp, E < 10⁻³ |
| E6 | Within-family trees | `phylo.py` | Poisson distance −ln(1−p), neighbour joining, rooted on a control |
| E7 | Motif search by multiple pattern matching | `motifs.py` | All PROSITE patterns that UniProt lists for the dataset, scanned at once as regexes |
| E8 | In-silico trypsin digest | `digest.py` | Cut after K/R unless the next residue is P; peptides of 6–30 aa; monoisotopic mass |
| E9 | Clustering and statistics | `analysis.py` | Standardized composition features; PCA; Ward clustering vs family (ARI); Kruskal–Wallis; sign test |
| E10 | AlphaFold confidence | `structure.py` | Mean pLDDT from AlphaFold DB |

---

## 5. Results

**E2/E4/E5: sequence similarity.** The headline result.

| Pair type | Pairs | Median identity | Median SW z-score |
|---|---|---|---|
| Same family | 194 | 68% | 136 |
| Luciferase vs its own control | 86 | 36% | 47 |
| Firefly vs click beetle (same ACS superfamily) | 32 | 49% | 126 |
| All other cross-family luciferase pairs | 1,259 | 20% | 0.0 (chance) |

- **BLAST.** There were 630 significant hits. The only cross-family hits link
  fireflies and click beetles (E ≤ 10⁻¹²⁶), which are expected relatives. No
  other luciferase found a luciferase from another family.
- **Weak alignment hits.** 36 of the 1,259 unrelated pairs reached z > 3
  (maximum 6.9). With 20 shuffles and over 1,000 tests, a few are expected
  by chance, and BLAST supported none of them.
- **Dot plots.** Same family: a continuous diagonal (479 dots). Luciferase vs
  control: a broken diagonal (334 dots). Firefly vs LuxA: noise only (6 dots).

**E3: dynamic programming.** The hand-written Needleman–Wunsch table gives a
score of −8, the same as Biopython's for identical settings
(`results/dp_example.txt`).

**E6: trees.** The trees reproduce known biology: firefly subfamilies,
*Vibrio* vs *Photobacterium* LuxA, and the *Armillaria* Luz cluster. In the
*Renilla*-type tree, the brittle-star luciferase groups with the **brittle
star's own non-luminous dehalogenase**, not with *Renilla* luciferase. This
directly illustrates two separate recruitments of the same enzyme family.

**E7: motifs.** The dataset carries 2 PROSITE patterns:
- **AMP-binding motif:** found in firefly and click-beetle luciferases **and**
  in their non-luminous controls. It is inherited from the family, not linked
  to light.
- **Bacterial-luciferase signature:** found in all 10 LuxA, and not in the
  SsuD controls.

No motif is shared between independent luciferase families.

**E8: peptides.** Of 1,806 tryptic peptides, only 2 occur in more than one
family, and both are shared by fireflies and click beetles. Each independent
family has its own peptide fingerprint, which could in principle identify it
by mass spectrometry.

**E1/E9: properties.**
- **Clustering.** Ward clustering of the composition features recovers the
  families with ARI = 0.80, where 1 is perfect. Most of the remaining mismatch
  is fireflies and click beetles mixing, as expected for relatives.
- **Family differences.** Every main property differs strongly between
  families (Kruskal–Wallis p < 10⁻⁵).

Luciferase vs own control, in the 5 families with controls:

| Property | Families where luciferase is… | Sign test p |
|---|---|---|
| GRAVY (hydrophobicity) | lower in 5 of 5 | 0.0625 |
| Aromaticity | higher in 5 of 5 | 0.0625 |
| Cysteine % | higher in 4 of 5 | 0.375 |
| pI | lower in 3, higher in 2 | 1.0 |

**E10: structure.** 62 of 67 proteins have AlphaFold models. The other five
are exactly the NCBI-only sequences. Median mean pLDDT is about 90 for the beetle, bacterial, fungal and *Renilla*
families, 86 for dinoflagellates and 56 for copepods. The low copepod value
suggests their small, disulfide-rich secreted luciferases are hard to
predict.

Figures are in `results/figures/`: `identity_heatmap.png`, `dotplots.png`,
`trees.png`, `pca_features.png` and `feature_boxplots.png`.

---

## 6. Conclusion

1. **H1 is supported: no shared sequence.** Independent luciferase families
   share no detectable sequence similarity, motifs or peptides. Alignment,
   dot plots, BLAST, PROSITE and peptide analysis agree. The only
   cross-family similarity (firefly with click beetle) reflects a shared
   ancestral enzyme family, not convergence.
2. **Properties are dominated by ancestry.** Composition sorts proteins by
   family (ARI 0.80), and the families differ strongly from each other.
3. **H2 has a consistent but unproven hint.** In all 5 testable families,
   luciferases are less hydrophobic and more aromatic than their non-luminous
   relatives. With 5 families each sign test gives p = 0.0625. It does not
   survive correction for testing 4 properties, and fireflies and click
   beetles are not fully independent. So this is a **hypothesis for future
   work**, not a finding.
4. **Convergence also happens within protein families.** Brittle star and
   *Renilla*, and fireflies and click beetles, each recruited the same enzyme
   family independently.

**Answer.** Independently evolved luciferases converged in **function**, not in
**sequence**. Their overall physicochemical profiles are set mostly by
ancestry, with a weak, consistent tendency towards lower hydrophobicity and
higher aromaticity that deserves testing on more families.

**Limitations.**
- Only 5 families have controls, and some groups are small (3 sequences).
- Several sequences are unreviewed or predicted, with evidence levels recorded.
- Secreted luciferases are analysed with their signal peptides.
- Whole-protein features cannot detect active-site or 3D-shape convergence.

**Future work.**
- Compare AlphaFold structures, for example the shape of the active-site
  pocket.
- Test the hydrophobicity and aromaticity hint on more luciferase families
  with controls.
- Analyse mature proteins, with signal peptides removed.

### Course coverage

| Week | Topic | Used in this project |
|---|---|---|
| 1 | Intro, origin of replication | Background only (genome-level topic; not applicable to protein comparison) |
| 2 | Motif search | PROSITE motif scan (E7) |
| 3 | Biological databases, GO | UniProt/NCBI curation, InterPro/Pfam, GO summary (`results/go_by_family.csv`) |
| 4 | Genome assembly | Not applicable (no genome data) |
| 5 | Antibiotics sequencing, mass spectrometry | Trypsin digest and peptide masses (E8) |
| 6 | Similarity and homology | Identity, z-scores and BLAST to separate homology from chance (E2, E5) |
| 7 | Global and local alignment | Needleman–Wunsch and Smith–Waterman (E2) |
| 8 | Dynamic programming | Hand-written DP table and traceback (E3) |
| 9 | Dot plots (sorting by reversals: not applicable) | Dot plots (E4) |
| 10 | Phylogenetics | Neighbour-joining trees (E6) |
| 11 | Clustering (transcriptomics: not applicable) | PCA, Ward clustering, ARI (E9) |
| 12 | Multiple pattern matching, BLAST | Multi-pattern PROSITE scan (E7), BLAST+ (E5) |
| 13 | Hidden Markov Models | Pfam/InterPro domains (profile-HMM based) used to confirm families; no HMM training |
| 14 | Proteomics, peptide sequencing | ProtParam properties (E1), tryptic peptides (E8) |
| 15 | Deep learning for biomolecules | AlphaFold structure confidence (E10) |

---

## References

1. Delroisse J, Duchatelet L, Flammang P, Mallefet J. Leaving the dark side? Insights into the evolution of luciferases. *Frontiers in Marine Science* 8:673620 (2021). https://doi.org/10.3389/fmars.2021.673620
2. He J, et al. Multiple origins of bioluminescence in beetles and evolution of luciferase function. *Molecular Biology and Evolution* 41(1):msad287 (2024). https://doi.org/10.1093/molbev/msad287
3. Kotlobay AA, et al. Genetically encodable bioluminescent system from fungi. *PNAS* 115(50):12728–12732 (2018). https://doi.org/10.1073/pnas.1803615115
4. Lau ES, et al. Functional characterization of luciferase in a brittle star indicates parallel evolution influenced by genomic availability of haloalkane dehalogenase. *Molecular Biology and Evolution* 42(5):msaf081 (2025). https://doi.org/10.1093/molbev/msaf081
5. Takenaka Y, et al. Evolution of bioluminescence in marine planktonic copepods. *Molecular Biology and Evolution* 29(6):1669–1681 (2012).
6. Liu L, Wilson T, Hastings JW. Molecular evolution of dinoflagellate luciferases, enzymes with three catalytic domains in a single polypeptide. *PNAS* 101(47):16555–16560 (2004).
7. Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology* 48(3):443–453 (1970).
8. Smith TF, Waterman MS. Identification of common molecular subsequences. *Journal of Molecular Biology* 147(1):195–197 (1981).
9. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. *PNAS* 89(22):10915–10919 (1992).
10. Camacho C, et al. BLAST+: architecture and applications. *BMC Bioinformatics* 10:421 (2009).
11. Saitou N, Nei M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4(4):406–425 (1987).
12. Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423 (2009).
13. Gasteiger E, et al. Protein identification and analysis tools on the ExPASy server. In: Walker JM (ed.), *The Proteomics Protocols Handbook*, Humana Press, 571–607 (2005).
14. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157(1):105–132 (1982).
15. Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589 (2021).
16. Varadi M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Research* 50(D1):D439–D444 (2022).
17. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51(D1):D523–D531 (2023).
18. Sigrist CJA, et al. New and continuing developments at PROSITE. *Nucleic Acids Research* 41(D1):D344–D347 (2013).
19. Kruskal WH, Wallis WA. Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association* 47(260):583–621 (1952).
20. Hubert L, Arabie P. Comparing partitions. *Journal of Classification* 2:193–218 (1985).

*Bibliographic details were written from memory, not checked against the
sources. Verify every reference before submission, especially [2]–[6].*

# Handoff Context: Bioinformatics Course Project

## 1. User and course context
- User is a Computer Science student taking Bioinformatics as an elective, not a biology specialist.
- This is a course project/presentation for the instructor, **not a plan to publish a paper**.
- The workload should be comparable to classmates' projects: not trivial, but not overly complex.
- The user prioritizes accessible public datasets, usable file formats, a coherent biological question, an implementable pipeline, interpretable results, and a clear conclusion.
- User wants explanations in simple language. Avoid overselling novelty or proposing unnecessarily large models, advanced wet-lab work, or complicated infrastructure.
- Proposal should fit in **no more than 2 pages in Microsoft Word**.

## 2. Chosen project
**Preferred full title:** “Molecular Signatures of Convergent Evolution: Comparative Sequence Analysis of Independently Evolved Luciferases”

**Short title alternative:** “Comparative Sequence Analysis of Independently Evolved Luciferases.”

**Main question:** Do independently evolved luciferases share common sequence-level or physicochemical characteristics despite having different evolutionary origins?

**Plain-language explanation:** We collect protein sequences of light-producing luciferases from fireflies, marine animals (e.g. *Renilla*), luminous bacteria, and fungi. We compare their sequences and measurable protein properties to see how similar or different they are, and whether different proteins that produce light show any common molecular characteristics.

**Important clarification:** This is a **protein-sequence comparative analysis** project. We are not primarily searching for a common DNA motif, and we are not training machine-learning models. We are implementing a Python bioinformatics pipeline, performing comparisons, visualizing results, and interpreting them.

## 3. Planned scope and method
Main luciferase groups:
1. Fireflies (firefly luciferase)
2. *Renilla* / sea pansies (Renilla luciferase)
3. Luminous bacteria (LuxA alpha subunit)
4. Bioluminescent fungi (Luz luciferase)

Where suitable sequences can be found, aim for around 5–10 representative proteins per group. An initial small starter dataset is useful for building and testing the pipeline, but group-level statistical conclusions are weak if only 1–3 sequences per group are used.

Pipeline:
1. Download annotated protein sequences from public resources such as UniProtKB and NCBI in **FASTA** format.
2. Curate the data: check protein names, organism, annotation, completeness, duplicates, and whether each sequence is actually a luciferase of the intended class.
3. Compute pairwise similarity. Potential tools: Biopython alignment functions and/or established alignment tools. Consider BLASTP because direct alignments between unrelated families can be difficult to interpret. Report identity/alignment scores carefully; low similarity alone does not prove independent evolutionary origin.
4. Extract protein features using Biopython: sequence length, molecular weight, amino acid composition, GRAVY/average hydrophobicity, aromaticity, and isoelectric point.
5. Make similarity matrices/heatmaps and plots/tables comparing feature distributions across families.
6. If suitable data are available, include selected non-luminous related proteins as controls to help distinguish inherited family properties from possible shared properties of light-producing proteins.
7. Use InterPro/domain annotations when useful to confirm protein-family backgrounds.

Suggested tools: Python, Biopython, Pandas, NumPy, Matplotlib, SciPy, UniProt, NCBI, InterPro. No GPU, deep learning, wet lab, or elaborate web app is required.

## 4. Expected deliverables
- Curated protein FASTA dataset, with accession numbers and source records documented.
- Python scripts/notebook that run the analysis reproducibly.
- Pairwise sequence similarity results and a heatmap/matrix.
- Protein physicochemical feature table and visualizations.
- Optional comparisons to relevant non-luminous proteins.
- Short interpretation, limitations, conclusions, and class presentation.

Expected patterns are hypotheses, not guaranteed results: proteins within a luciferase family may be more similar to one another than to proteins from other families. The cross-family physicochemical analysis tests whether any common patterns exist; do not assume or claim that such patterns must exist or cause bioluminescence.

## 5. Important biology and scientific caveats
- Bioluminescence evolved multiple times. Many luciferase systems come from different protein families, but do not claim every luciferase is unrelated to every other luciferase.
- Firefly luciferases are related to the acyl-CoA synthetase superfamily.
- *Renilla* luciferase and brittle-star luciferase involve independent recruitment from proteins related to the haloalkane dehalogenase family; this is an example of parallel evolution from a related ancestral family, so these should not be described as wholly unrelated.
- Fungal Luz is a distinct luciferase system.
- Pairwise percent identity by itself cannot prove an absence of homology. Avoid a single phylogenetic tree that forces all very different luciferase families into one alignment/tree. Use family annotations and published evolutionary evidence to contextualize results.
- Comparing a luminous protein with one non-luminous protein is suggestive, not enough to prove that a given property is caused by light-producing function. Use cautious conclusions.
- Genomic clustering, DNA GC content, codon usage, and neural models are not part of the agreed core scope; only add them if there is a clear need and time.

## 6. Candidate starter sequences mentioned in prior discussion
These are **candidate accessions from earlier planning and must be rechecked against UniProt before being treated as a finalized dataset**. Verify current record, organism, protein name, annotation status, and whether sequence is complete.

Firefly candidates:
- *Photinus pyralis*, P08659: https://rest.uniprot.org/uniprotkb/P08659.fasta
- *Nipponoluciola cruciata*, P13129: https://rest.uniprot.org/uniprotkb/P13129.fasta
- *Lampyris turkestanicus*, Q5UFR2: https://rest.uniprot.org/uniprotkb/Q5UFR2.fasta

Renilla candidates:
- *Renilla reniformis*, P27652: https://rest.uniprot.org/uniprotkb/P27652.fasta
- *Renilla muelleri*, Q9BLZ3: https://rest.uniprot.org/uniprotkb/Q9BLZ3.fasta

Bacterial LuxA candidates:
- *Aliivibrio fischeri*, P19907: https://rest.uniprot.org/uniprotkb/P19907.fasta
- *Photobacterium leiognathi*, P09140: https://rest.uniprot.org/uniprotkb/P09140.fasta

Fungal Luz candidates:
- *Neonothopanus nambi*, A0A3G9JYH7: https://rest.uniprot.org/uniprotkb/A0A3G9JYH7.fasta
- *Armillaria gallica*, A0A2H3E985: https://rest.uniprot.org/uniprotkb/A0A2H3E985.fasta

Possible non-luminous comparison proteins:
- *Drosophila melanogaster* fatty acyl-CoA synthetase CG6178, Q9VCC6: https://rest.uniprot.org/uniprotkb/Q9VCC6.fasta
- *Xanthobacter autotrophicus* haloalkane dehalogenase DhlA, P22643: https://rest.uniprot.org/uniprotkb/P22643.fasta

Notes: The controls above are starting candidates, not automatically the closest biological homologs. In particular, DhlA is a representative haloalkane dehalogenase, not necessarily the closest non-luminous homolog of *Renilla*. Find/verify suitable homologs before making strong evolutionary comparisons. Do not assume the accession list or a combined download link is correct without checking records.

UniProt FASTA REST pattern: `https://rest.uniprot.org/uniprotkb/ACCESSION.fasta`
UniProt page pattern: `https://www.uniprot.org/uniprotkb/ACCESSION`
NCBI Protein: https://www.ncbi.nlm.nih.gov/protein/
InterPro: https://www.ebi.ac.uk/interpro/

## 7. Core literature / references previously used
1. J. Delroisse, L. Duchatelet, P. Flammang, and J. Mallefet, “Leaving the Dark Side? Insights Into the Evolution of Luciferases,” *Frontiers in Marine Science*, vol. 8, Art. no. 673620, 2021. https://doi.org/10.3389/fmars.2021.673620
   - Relevant broad comparison of luciferases/photoproteins; previous discussion noted it covers 134 sequences from 75 species across 11 phyla. This means the broad idea of comparing luciferase evolutionary groups is not wholly new, which is fine for this course project.
2. J. He et al., “Multiple Origins of Bioluminescence in Beetles and Evolution of Luciferase Function,” *Molecular Biology and Evolution*, vol. 41, no. 1, Art. no. msad287, 2024. https://doi.org/10.1093/molbev/msad287
3. A. A. Kotlobay et al., “Genetically Encodable Bioluminescent System from Fungi,” *Proceedings of the National Academy of Sciences*, vol. 115, no. 50, pp. 12728–12732, 2018. https://doi.org/10.1073/pnas.1803615115
4. E. S. Lau et al., “Functional Characterization of Luciferase in a Brittle Star Indicates Parallel Evolution Influenced by Genomic Availability of Haloalkane Dehalogenase,” *Molecular Biology and Evolution*, vol. 42, no. 5, Art. no. msaf081, 2025. https://doi.org/10.1093/molbev/msaf081

## 8. Proposal status and approved structure
The user has a compact proposal intended to fit within 2 Word pages, with these headings in this order:
1. Abstract
2. Introduction
3. Literature Review
4. Proposed Methodology / Tentative Algorithms and Tools
5. Expected Outcome or Deliverables
6. References

Proposal content is broadly approved, but use the following small accuracy improvements if editing it:
- Say “Many luciferase systems evolved independently…” rather than claiming all different luciferases necessarily evolved independently from unrelated families.
- Methodology can mention pairwise alignment using Needleman–Wunsch/Smith–Waterman and/or BLASTP; BLASTP helps assess similarity without forcing an alignment between very different proteins.
- Keep the conclusion cautious: “The analysis will determine whether independently evolved luciferases remain distinct at the sequence level while showing any common physicochemical patterns despite their different evolutionary origins and shared light-producing function.” Avoid implying the analysis can prove a physicochemical property causes light production.
- The abstract is around 170 words in the shortest proposal version, which is okay if 250 words is a maximum, but would need expansion if the instructor means approximately/exactly 250.

## 9. Short project description (ready to reuse)
“This project investigates whether independently evolved luciferase proteins from fireflies, marine animals, bacteria, and fungi share common molecular characteristics despite originating from different protein families. Using publicly available protein sequences from UniProt and NCBI, we will perform sequence alignment, similarity analysis, and physicochemical feature extraction such as protein length, amino acid composition, hydrophobicity, molecular weight, aromaticity, and isoelectric point. Selected luciferases will also be compared with related non-luminous proteins, and the results will be visualized using heatmaps, graphs, and comparison tables to study possible patterns of convergent evolution.”

## 10. The user’s latest concise explanation
“We’re comparing luciferase protein sequences and properties across fireflies, bacteria, fungi, and marine animals to see whether independently evolved proteins share common molecular characteristics.”

## 11. How to help next
- Keep recommendations focused and implementable for an undergraduate CS elective.
- If continuing data work, first validate each candidate UniProt accession from the live records and gather a final balanced dataset (or at least make it clear that the starter set is small).
- Give direct, clickable dataset links and explain file formats simply.
- The best immediate next step is likely verifying the proposed FASTA accessions and preparing a clean data manifest: family, species, protein name, accession, URL, annotation status, sequence length, and why included.

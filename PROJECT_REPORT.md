# Project Report: Comparing Independently Evolved Luciferases

A plain-language record of what we are doing, why, how, and what we found.
It is updated after every phase. Last updated: 2026-11-04.

| Phase | Status |
|---|---|
| 1. Build the dataset | Done |
| 2. Set up BLAST | Done |
| 3. Calculate protein properties | Done |
| 4. Alignment and similarity heatmap | Not started |
| 5. Dot plots | Not started |
| 6. BLAST searches | Not started |
| 7. Family trees | Not started |
| 8. Clustering and control comparison | Not started |
| 9. Notebook, conclusions, presentation | Not started |

---

## The project in brief

**Background.** Many living things make their own light: fireflies, some
bacteria, mushrooms, tiny sea creatures. The protein that drives the
light-making reaction is called a **luciferase**. Light production evolved
many separate times, so these luciferases are often not related to each
other, even though they do a similar job.

**Our question.** Do independently evolved luciferases share any common
features in their protein sequences or physical and chemical properties,
despite their different origins?

**Why it is interesting.** When unrelated organisms arrive at a similar
solution, it is called **convergent evolution**. We want to see whether that
leaves a measurable mark on the proteins.

**What this project is not.**
- It compares **proteins**, not DNA.
- It is computer analysis of public data, with no lab work.
- It does not use machine learning.
- It cannot prove that a property *causes* light production. It can only
  show patterns.

---

## Phase 1: Build the dataset (done)

**Purpose.** Collect a reliable set of protein sequences to analyse.

**Why it matters.** Every later result depends on this data. A wrongly
labelled or incomplete protein would give misleading results.

**What we did.**
1. Searched the public databases UniProt and NCBI for luciferases.
2. Compared candidates with each other to remove near-copies.
3. Added **controls** to the first four families.
4. Downloaded everything with a script, then checked it with a second script.

**What is a control?** A non-light-producing "cousin" of a luciferase, from
the same protein family. If a luciferase has some property, the control shows
whether that property is linked to making light or is simply normal for that
family. If both have it, it came from the family. If only the luciferase has
it, it might be linked to light.

**Rules for including a sequence.**
- There must be good reason to believe it really is a luciferase.
- It must be complete. Six sequences that may lack a few end residues were
  kept and flagged.
- Sequences within a family must be under 97% identical, so that one group of
  close relatives does not dominate.
- Expert-reviewed database entries are preferred.

**Result: 58 proteins (51 luciferases, 7 controls).**

| Family | Luciferases | Controls | Example species |
|---|---|---|---|
| Firefly | 8 | 2 | *Photinus*, *Luciola*, *Photuris* |
| *Renilla*-type | 3 | 3 | *Renilla* (sea pansy), brittle star |
| Bacterial LuxA | 10 | 1 | *Vibrio*, *Photobacterium* |
| Fungal Luz | 9 | 1 | *Neonothopanus*, *Armillaria*, *Mycena* |
| Copepod | 10 | 0 | *Gaussia*, *Metridia* |
| Dinoflagellate | 8 | 0 | *Lingulodinium*, *Pyrocystis* |
| Ostracod | 3 | 0 | *Vargula*, *Cypridina* |

**Why seven families?** We started with four. The main question is about
independent origins, so each new family adds more evidence than extra
sequences within a family we already have. The sea pansy and the brittle star
evolved their luciferases separately, so the seven families hold eight origins.

**Limits of this phase.**
- *Renilla*-type and ostracod have only 3 sequences each, because that is all
  the public databases offer with good evidence.
- The copepod, dinoflagellate and ostracod families have no controls, because
  no suitable non-luminous relatives are known.
- Only 21 of the 58 entries are expert-reviewed. The evidence level of each
  entry is recorded in the manifest.
- That the controls produce no light of their own comes from the literature.
  We still need to cite sources for this.

**Where to find it.** `data/curated/` holds the sequences and
`manifest.csv`. `data/README.md` has the full rules.

---

## Phase 2: Set up BLAST (done)

**Purpose.** Install BLAST, the standard tool for finding similar sequences.

**Why.** NCBI's website version was too slow for repeated use; one search was
still waiting after 10 minutes. Running it on our own computer takes seconds
and is easy to repeat.

**What we did.** Installed BLAST+ 2.17.0 and ran one test search, using a
firefly luciferase against all 58 sequences.

**Result of the test.** It matched only the other 7 fireflies (60–85%
identical) and the 2 firefly controls (42–46%). It found nothing in the other
six families. This is only one search, but it is the pattern we expected.

---

## Phase 3: Calculate protein properties (done)

**Purpose.** Measure each protein so that families can be compared with
numbers, even when their sequences look nothing alike.

**Why.** Unrelated proteins cannot be compared letter by letter, but all of
them have a size, a charge and a hydrophobicity. This is how we test for
shared features.

**What we measured.**

| Property | Plain meaning |
|---|---|
| Length | Number of amino acids in the chain |
| Molecular weight | How heavy the protein is, in kilodaltons (kDa) |
| Isoelectric point (pI) | The pH at which the protein has no overall charge |
| GRAVY | Average hydrophobicity; more negative means more water-loving |
| Aromaticity | Share of ring-shaped amino acids (F, W, Y) |
| Amino acid composition | Percentage of each of the 20 amino acids |

**How.** The script `src/features.py` uses the ProtParam tool in the Biopython
library. To check it, we confirmed there were no missing values, that the
amino acid percentages add up to 100, and that two well-known proteins came out
at sensible weights.

**Results (median of each group).**

| Family | Group | n | Length | Weight (kDa) | pI | GRAVY | Aromaticity |
|---|---|---|---|---|---|---|---|
| Firefly | luciferase | 8 | 547.5 | 60.5 | 6.45 | -0.080 | 0.093 |
| | control | 2 | 544 | 60.0 | 8.42 | -0.016 | 0.081 |
| *Renilla*-type | luciferase | 3 | 311 | 36.1 | 5.86 | -0.430 | 0.119 |
| | control | 3 | 296 | 33.2 | 5.03 | -0.228 | 0.109 |
| Bacterial LuxA | luciferase | 10 | 354.5 | 40.4 | 5.38 | -0.378 | 0.104 |
| | control | 1 | 381 | 41.7 | 5.69 | -0.192 | 0.079 |
| Fungal Luz | luciferase | 9 | 261 | 29.7 | 8.63 | -0.236 | 0.086 |
| | control | 1 | 283 | 31.4 | 8.16 | -0.172 | 0.092 |
| Copepod | luciferase | 10 | 210.5 | 23.0 | 7.08 | -0.101 | 0.053 |
| Dinoflagellate | luciferase | 8 | 1,237 | 136.9 | 6.16 | -0.491 | 0.101 |
| Ostracod | luciferase | 3 | 553 | 61.5 | 4.66 | -0.211 | 0.097 |

**Observations (early, not conclusions).**
1. **No single "luciferase profile" so far.** Size ranges from about 210 to
   1,237 amino acids, and pI from 4.66 to 8.63.
2. **Hydrophobicity is worth a closer look.** In all 4 families with controls,
   the luciferase is less hydrophobic than its control. This could still be
   chance: four out of four has a probability of about 0.125, some controls
   are a single protein, and we looked at several properties.
3. **Aromaticity** is higher than the control in 3 of 4 families.
4. **pI shows no pattern.** The difference from the control goes in both
   directions.
5. **Copepod and ostracod luciferases are rich in cysteine** (about 5–6%,
   against 1–2.5% elsewhere). Both are secreted proteins, so this may be
   related to secretion rather than to making light.

**Limits of this phase.**
- Length and weight mostly reflect protein structure. Dinoflagellate
  luciferases are large because they contain three repeated domains. GRAVY,
  pI, aromaticity and composition are fairer comparisons.
- Secreted luciferases still include their signal peptide, a short tag removed
  from the mature protein. This shifts their values slightly.
- Proteins in one family share ancestry, so they are not independent samples.
  The family, not the single protein, is the real unit of evidence.

**Where to find it.** `results/features.csv` has one row per protein.
`results/features_by_family.csv` has family summaries. `results/README.md`
explains each column.

---

## What comes next

- **Alignment and heatmap:** measure how similar every pair of proteins is. We
  expect high similarity within families and very little between them.
- **Dot plots:** show similarity between two proteins as a picture.
- **BLAST searches:** confirm that luciferases find only their own relatives.
- **Family trees:** show relationships within each family. We will not build
  one tree of all families, because unrelated proteins cannot be fairly
  aligned together.
- **Clustering and control comparison:** check whether proteins group by family
  or by light production, and test the hydrophobicity pattern properly.

---

## Glossary

| Term | Meaning |
|---|---|
| Amino acid / residue | One building block of a protein; a protein is a chain of them |
| Luciferase | A protein that drives a light-producing reaction |
| Convergent evolution | Unrelated organisms independently evolving a similar trait |
| Family | A group of proteins that share a common ancestor |
| Homolog | A protein related to another by shared ancestry |
| Control | A related protein that does not make light, used for comparison |
| FASTA | A simple text format for storing sequences |
| UniProt / NCBI | Public databases of biological sequences |
| Alignment | Lining up two sequences to see where they match |
| BLAST | A tool that finds sequences similar to a query |
| Median | The middle value; less affected by extreme values than the average |

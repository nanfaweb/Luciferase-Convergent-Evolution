# Dataset

37 protein sequences: 30 luciferases from four independently evolved systems,
plus 7 non-luminous homologs used as controls. Sources are UniProtKB (34) and
NCBI Protein (3). Retrieved 2026-11-03.

| Family | Luciferases | Controls (same protein family, non-luminous) | Shared Pfam domain |
|---|---|---|---|
| `firefly` | 8 (5 lampyrid subfamilies) | AbLL from a non-luminous click beetle (A9ZPM4), *Drosophila* CG6178 (Q9VCC6) | PF00501 AMP-binding |
| `renilla_type` | 3: *R. reniformis*, *R. muelleri*, brittle star *Amphiura filiformis* (AfLuc) | Bacterial dehalogenases LinB (D4Z2G1), DhaA (P0A3G2); brittle-star dehalogenase DafA (XRP75329.1) | PF00561 Abhydrolase_1 |
| `bacterial_luxa` | 10 (*Aliivibrio*, *Vibrio*, *Photobacterium*, *Photorhabdus*, flashlight-fish symbionts) | Alkanesulfonate monooxygenase SsuD (P80645) | PF00296 Bac_luciferase |
| `fungal_luz` | 9 (*Neonothopanus* ×2, *Armillaria* ×2, *Mycena* ×3, *Omphalotus*, *Panellus*) | Luz-domain protein from non-luminous *Fistulina hepatica* (A0A0D7AKP7) | PF17648 Luciferase |

In each family, the UniProt controls share the Pfam domain of their
luciferases. The domain annotations come from UniProt and confirm that each
family assignment is correct. NCBI records carry no Pfam cross-references.

`renilla_type` covers luciferases recruited from the haloalkane dehalogenase
(HLD) family. *Renilla* and the brittle star recruited it **independently**:
AfLuc is only 45% identical to *Renilla* luciferase. So this group holds two
lineages, *Renilla* (2) and *Amphiura* (1). That parallel recruitment is
itself a useful within-group comparison.

## Files

| Path | Format | What it is |
|---|---|---|
| `seeds.csv` | CSV | **Source of truth.** One row per accession: `source` (`uniprot`/`ncbi`), `family`, `role` (`luciferase`/`control`), optional `label` and `organism` overrides, `allow_partial`, `strip_suffix`, and the `reason` it was included. Edit this file to change the dataset. |
| `raw/{ACCESSION}.fasta` | FASTA | The record exactly as the database served it (tags included). |
| `curated/manifest.csv` | CSV | One row per sequence: source, organism, taxon ID, reviewed, evidence level, curated length, fragment / partial-accepted / residues-stripped / complete flags, entry and sequence versions, SHA-256 of the curated sequence, function note, GO molecular function, InterPro, Pfam, URLs, retrieval date. |
| `curated/luciferases.fasta` | FASTA | The 30 curated luciferases. Header: `>ACC|family|role|organism|label` |
| `curated/controls.fasta` | FASTA | The 7 curated controls, same header format. |
| `curated/all_sequences.fasta` | FASTA | All 37, in `seeds.csv` order. **Use the curated files for analysis, not `raw/`.** |
| `curated/validation_summary.json` | JSON | Counts per family and role, completeness, reviewed count. |

Rebuild and check:

```bash
python -m src.download
python -m src.verify_data --live
```

`download` fetches and checks every record before it writes anything. It
also deletes raw FASTA files that are no longer listed in `seeds.csv`.
`verify_data` checks the files against the seeds, sequence hashes, residues,
the combined FASTA files and the redundancy rule. With `--live`, it also
checks UniProt records for updates. NCBI accessions are versioned, so a
changed sequence would get a new accession.

## Selection rules

1. **Functional evidence.** A sequence counts as a luciferase only if it is a
   characterised luciferase, or an ortholog from a known luminous species
   within ortholog-level identity of characterised ones (≥65% for Luz, ≥77%
   for LuxA). Many genome-annotated "alkanal monooxygenase" and "luciferase
   domain" proteins in non-luminous bacteria and fungi were rejected. So were
   genome-annotated "Renilla-type luciferases" (e.g. *Galeolaria*,
   *Paramuricea*) with no functional data.
2. **Complete sequences.** The sequence must start with Met. Entries flagged
   as fragments were rejected, with one documented exception. The five Luz
   cDNA clones from Kotlobay et al. 2018 (*Omphalotus*, *Panellus*, *N.
   gardneri*, *M. citricolor*, *M. chlorophos*) have `allow_partial=yes`.
   UniProt marks only their last residue as non-terminal, because the clones
   stop before the stop codon. Compared with full-length *N. nambi* Luz, they
   lack 0–6 C-terminal residues. Their length and molecular weight may
   therefore be slightly underestimated; the `partial_accepted` column marks
   them.
3. **Expression tags.** AfLuc and DafA come from the GenBank expression
   constructs deposited by Lau et al. 2025 (nucleotide PP777633 and
   PP777634). The C-terminal His6 tag is removed (`strip_suffix`), and the
   organism is set to *A. filiformis* because GenBank lists "synthetic
   construct". No native RefSeq protein for either gene could be found.
4. **Evidence.** Reviewed (Swiss-Prot) entries are preferred. Unreviewed
   entries add genus or subfamily breadth, and `evidence_level` records how
   strong each one is.
5. **Redundancy.** At most one sequence per species, except *P. leiognathi*,
   which has two distinct LuxA genes (77% identical). No two luciferases in
   the same family may be ≥97% identical (global alignment, BLOSUM62).
   Dropped by this rule: *Photorhabdus* C7BN05 (98% to P23146), *P.
   kishitanii* Q6VFP7 (100% to P24113), *A. salmonicida* B6ESM5 (98% to
   *A. logei*), and the *M. chlorophos* genome prediction A0A8H6SUE6
   (replaced by its 98.8%-identical transcript clone). Also left out, to
   avoid over-sampling one genus while staying within 10 per family:
   *Pyrocoelia* Q9GPF9 and A1Z0H9 (95–97% to *Lampyris*), *Photorhabdus*
   P19839 and Q7N575 (93% to P23146), *V. chagasii* B2XS59 (94% to *V.
   harveyi*), and *Armillaria ostoyae* and *A. fuscipes* (95–96% to *A.
   gallica*).
6. **Controls.** Each control is the closest non-luminous homolog found by
   global identity to its family's luciferases (see `reason` in `seeds.csv`).

## Known limitations (state these in the write-up)

- **`renilla_type` has only 3 luciferases from 2 lineages.** This is the
  limit of public data, not of the search. UniProt and NCBI have no other
  functionally characterised HLD-derived luciferases; *Pyrosoma* PyroLuc was
  left out because Lau et al. report only weak activity. Report this group
  descriptively, with no group-level statistics.
- **Unequal sample sizes:** 8, 3, 10 and 9 luciferases. Use per-family
  medians and non-parametric tests, and do not pool everything into one test.
- **Phylogenetic non-independence:** sequences within a family share
  ancestry, so they are not independent samples. Conclusions about
  "luciferase vs control" should rest on the 4 family-level comparisons, not
  on 30 sequences.
- **Unchecked literature claims.** That AbLL (A9ZPM4), CG6178 (Q9VCC6) and
  brittle-star DafA give off no native light comes from the literature and
  was not re-checked here. UniProt names Q9VCC6 and A9ZPM4 "Luciferin
  4-monooxygenase" only by homology. Cite sources before relying on these
  claims. The *Fistulina* control is a genome prediction.
- **Evidence levels for LuxA.** P19907, P09140 and P18299 are "inferred from
  homology" at the protein-existence level, even though they are reviewed
  entries. The three unreviewed LuxA additions are also homology-inferred;
  they are luxA genes from luminous species.
- Sequence similarity alone cannot rule out deep homology across families.
  See the caveats in `../bioinformatics_project_handoff.md` §5.

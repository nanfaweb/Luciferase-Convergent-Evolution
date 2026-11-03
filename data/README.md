# Dataset

27 UniProtKB protein sequences: 21 luciferases from four independently evolved
systems, plus 6 non-luminous homologs used as controls. Retrieved 2026-11-03.

| Family | Luciferases | Controls (same protein family, non-luminous) | Shared Pfam domain |
|---|---|---|---|
| `firefly` | 8 (5 lampyrid subfamilies) | AbLL from a non-luminous click beetle (A9ZPM4), *Drosophila* CG6178 (Q9VCC6) | PF00501 AMP-binding |
| `renilla` | 2 (*R. reniformis*, *R. muelleri*) | Haloalkane dehalogenases LinB (D4Z2G1), DhaA (P0A3G2) | PF00561 Abhydrolase_1 |
| `bacterial_luxa` | 7 (Vibrionaceae, *Photorhabdus*, fish symbiont) | Alkanesulfonate monooxygenase SsuD (P80645) | PF00296 Bac_luciferase |
| `fungal_luz` | 4 (*Neonothopanus*, 2× *Armillaria*, *Mycena*) | Luz-domain protein from non-luminous *Fistulina hepatica* (A0A0D7AKP7) | PF17648 Luciferase |

In each family, the controls share the Pfam domain of their luciferases. The
domain annotations come from UniProt and confirm that each family assignment
is correct.

## Files

| Path | Format | What it is |
|---|---|---|
| `seeds.csv` | CSV | **Source of truth.** One row per accession: `family`, `role` (`luciferase`/`control`), optional display `label`, and the `reason` it was included. Edit this file to change the dataset. |
| `raw/{ACCESSION}.fasta` | FASTA | The record exactly as UniProt served it (`>sp|ACC|NAME ...` header). |
| `curated/manifest.csv` | CSV | One row per sequence with UniProt metadata: organism, taxon ID, reviewed, evidence level, length, fragment/complete flags, entry and sequence versions, SHA-256 of the sequence, function note, GO molecular function, InterPro, Pfam, URLs, retrieval date. |
| `curated/luciferases.fasta` | FASTA | The 21 luciferases. Header: `>ACC|family|role|organism|label` |
| `curated/controls.fasta` | FASTA | The 6 controls, same header format. |
| `curated/all_sequences.fasta` | FASTA | All 27, in `seeds.csv` order. |
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
checks for UniProt updates since the download.

## Selection rules

1. **Complete sequences only.** Entries marked as UniProt fragments, entries
   with non-terminal residues, and sequences that do not start with Met are
   rejected. This excluded the transcript-level Luz sequences from
   *Panellus*, *Omphalotus* and *Mycena citricolor*, which are all fragments.
2. **Evidence.** Reviewed (Swiss-Prot) entries are preferred. Unreviewed
   entries are used only when they add a new genus or subfamily, and the
   `evidence_level` column records how strong each one is.
3. **Redundancy.** At most one sequence per species, except
   *P. leiognathi*, which has two distinct LuxA genes (77% identical). No two
   luciferases in the same family may be ≥97% identical (global alignment,
   BLOSUM62). Dropped for this reason: *Pyrocoelia* Q9GPF9 and A1Z0H9 (≥95%
   to *Lampyris*), *Photorhabdus* P19839 and Q7N575 (≥93% to P23146),
   *Armillaria ostoyae* (96% to *A. gallica*).
4. **Controls.** Each control is the closest non-luminous homolog found by
   global identity to its family's luciferases (see `reason` in `seeds.csv`).
   The controls replace the earlier choice, *Xanthobacter* DhlA (P22643). DhlA
   is only 25% identical to *Renilla* luciferase, while LinB and DhaA are
   36–40% identical.

## Known limitations (state these in the write-up)

- **The *Renilla* group has two sequences, and they are 96.5% identical.**
  UniProt and NCBI have no other functionally characterised Renilla-type
  luciferases. The brittle-star luciferase (Lau et al., 2025) is not
  deposited under a searchable protein name. Treat this group as a single
  lineage, and do not report group-level statistics for it.
- **Fungal Luz** has 4 sequences, and two of them are from *Armillaria*.
  The *Mycena chlorophos* sequence is a genome prediction (evidence level 4).
- **Labels for controls.** UniProt names Q9VCC6 and A9ZPM4 "Luciferin
  4-monooxygenase", but that name comes only from automatic homology
  annotation. `seeds.csv` gives them functional labels instead. Studies of
  these beetle and fly homologs describe them as acyl-CoA synthetases with no
  native light emission. [This claim comes from the literature and was not
  re-checked here; cite a source before using it in the report.] For the same
  reason, the *Fistulina* Luz-domain control is only predicted.
- **Evidence levels for LuxA.** P19907, P09140 and P18299 are "inferred
  from homology" at the protein-existence level, even though they are
  reviewed entries.
- Sequence similarity alone cannot rule out deep homology across families.
  See the caveats in `../bioinformatics_project_handoff.md` §5.

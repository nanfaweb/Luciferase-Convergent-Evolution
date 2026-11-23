"""In-silico trypsin digest (W5 mass spectrometry, W14 peptide sequencing).

Trypsin cuts after K or R, except when the next residue is P. Peptides of 6-30
residues (the range a mass spectrometer usually detects) get a monoisotopic
mass. If independently evolved luciferases shared sequence, they would share
peptides; peptides unique to one family could identify it by mass.

Writes results/tryptic_peptides.csv and results/tryptic_summary.csv.

Usage:
    python -m src.digest
"""

from __future__ import annotations

import re

import pandas as pd
from Bio.SeqUtils import molecular_weight

from src.align import load_proteins
from src.features import RESULTS_DIR

MIN_LEN, MAX_LEN = 6, 30


def trypsin_digest(sequence: str) -> list[str]:
    return [p for p in re.split(r"(?<=[KR])(?!P)", sequence) if p]


def run() -> pd.DataFrame:
    proteins = load_proteins()
    rows = []
    for prot in proteins.itertuples():
        for peptide in trypsin_digest(prot.sequence):
            if MIN_LEN <= len(peptide) <= MAX_LEN:
                rows.append({"accession": prot.accession, "family": str(prot.family), "role": prot.role,
                             "peptide": peptide,
                             "monoisotopic_mass": molecular_weight(peptide, "protein", monoisotopic=True)})
    peptides = pd.DataFrame(rows)
    families_per_peptide = peptides.groupby("peptide")["family"].nunique()
    peptides["n_families_with_peptide"] = peptides["peptide"].map(families_per_peptide)
    peptides.to_csv(RESULTS_DIR / "tryptic_peptides.csv", index=False, float_format="%.4f", lineterminator="\n")

    summary = peptides.groupby(["family", "role"], sort=False).agg(
        proteins=("accession", "nunique"), peptides=("peptide", "size"),
        shared_with_other_family=("n_families_with_peptide", lambda n: int((n > 1).sum())),
    ).reset_index()
    summary.to_csv(RESULTS_DIR / "tryptic_summary.csv", index=False, lineterminator="\n")
    shared = peptides[peptides["n_families_with_peptide"] > 1]
    print(f"{len(peptides)} tryptic peptides ({MIN_LEN}-{MAX_LEN} aa); "
          f"{shared['peptide'].nunique()} distinct peptides occur in more than one family")
    if len(shared):
        print(shared.groupby("peptide")["family"].agg(lambda f: ", ".join(sorted(set(f)))).to_string())
    return summary


if __name__ == "__main__":
    run()

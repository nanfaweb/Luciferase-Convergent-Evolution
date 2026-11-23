"""AlphaFold (deep-learning) structure confidence per protein (W15).

Queries the AlphaFold Protein Structure Database for each UniProt accession and
records the mean pLDDT (0-100; above 70 means the predicted fold is generally
reliable, above 90 very reliable). NCBI-only sequences have no AlphaFold entry.
This shows how well predicted 3D structures cover the dataset, the starting
point for comparing luciferases at the structure level.

Writes results/alphafold_confidence.csv.

Usage:
    python -m src.structure
"""

from __future__ import annotations

import time

import pandas as pd
import requests

from src.align import load_proteins
from src.features import RESULTS_DIR

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def run() -> pd.DataFrame:
    proteins = load_proteins()
    session = requests.Session()
    rows = []
    for prot in proteins.itertuples():
        row = {"accession": prot.accession, "family": str(prot.family), "role": prot.role,
               "alphafold_id": "", "mean_plddt": None, "sequence_matches": None}
        if prot.source == "uniprot":
            response = session.get(AFDB_API.format(accession=prot.accession), timeout=60)
            if response.ok and response.json():
                model = response.json()[0]
                row.update(alphafold_id=model.get("entryId", ""),
                           mean_plddt=model.get("globalMetricValue"),
                           sequence_matches=model.get("uniprotSequence") == prot.sequence)
            time.sleep(0.2)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS_DIR / "alphafold_confidence.csv", index=False, float_format="%.1f", lineterminator="\n")
    found = table.dropna(subset=["mean_plddt"])
    print(f"AlphaFold models for {len(found)} of {len(table)} proteins")
    print(found.groupby(["family", "role"], sort=False)["mean_plddt"].agg(["count", "median", "min"]).round(1))
    return table


if __name__ == "__main__":
    run()

"""Calculate physicochemical properties for every curated protein.

Reads data/curated/all_sequences.fasta and manifest.csv, then writes:

- results/features.csv: one row per protein
- results/features_by_family.csv: n, median, min and max per family and role

Usage:
    python -m src.features
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from src.download import CURATED_DIR, FAMILIES, PROJECT_ROOT
from src.verify_data import parse_fasta

RESULTS_DIR = PROJECT_ROOT / "results"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ID_COLUMNS = ["accession", "family", "role", "organism", "label"]
FLAG_COLUMNS = ["reviewed", "partial_accepted"]
SUMMARY_FEATURES = ["length", "molecular_weight_kda", "isoelectric_point", "gravy", "aromaticity"]


def protein_features(sequence: str) -> dict:
    analysis = ProteinAnalysis(sequence)
    features = {
        "length": len(sequence),
        "molecular_weight_kda": analysis.molecular_weight() / 1000,
        "isoelectric_point": analysis.isoelectric_point(),
        "gravy": analysis.gravy(),
        "aromaticity": analysis.aromaticity(),
    }
    composition = analysis.amino_acids_percent
    for residue in AMINO_ACIDS:
        features[f"aa_{residue}_pct"] = composition[residue]
    return features


def build_feature_table(fasta_path: Path, manifest_path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, keep_default_na=False)
    records = parse_fasta(fasta_path.read_text(encoding="utf-8"))
    features = pd.DataFrame(
        [{"accession": acc, **protein_features(seq)} for acc, seq in records]
    )
    table = manifest[ID_COLUMNS + FLAG_COLUMNS].merge(
        features, on="accession", how="outer", validate="one_to_one", indicator=True
    )
    if not (table.pop("_merge") == "both").all():
        raise ValueError("manifest and FASTA list different accessions")
    order = {acc: i for i, acc in enumerate(manifest["accession"])}
    return table.sort_values("accession", key=lambda s: s.map(order)).reset_index(drop=True)


def summarise(table: pd.DataFrame) -> pd.DataFrame:
    grouped = table.groupby(["family", "role"], sort=False)
    summary = grouped[SUMMARY_FEATURES].agg(["median", "min", "max"])
    summary.columns = [f"{feature}_{stat}" for feature, stat in summary.columns]
    summary.insert(0, "n", grouped.size())
    summary = summary.reset_index()
    summary["family"] = pd.Categorical(summary["family"], FAMILIES, ordered=True)
    return summary.sort_values(["family", "role"], ascending=[True, False])


def run() -> pd.DataFrame:
    RESULTS_DIR.mkdir(exist_ok=True)
    table = build_feature_table(
        CURATED_DIR / "all_sequences.fasta", CURATED_DIR / "manifest.csv"
    )
    table.to_csv(
        RESULTS_DIR / "features.csv", index=False, float_format="%.4f", lineterminator="\n"
    )
    summarise(table).to_csv(
        RESULTS_DIR / "features_by_family.csv",
        index=False,
        float_format="%.3f",
        lineterminator="\n",
    )
    print(f"wrote {len(table)} proteins to results/features.csv")
    return table


if __name__ == "__main__":
    run()

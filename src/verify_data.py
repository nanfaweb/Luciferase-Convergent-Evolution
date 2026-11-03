"""Verify the curated dataset against data/seeds.csv, the local FASTA files and
(with --live) the current UniProt records.

Usage:
    python -m src.verify_data          # offline checks
    python -m src.verify_data --live   # also compare with live UniProt
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import time

import requests
from Bio.Align import PairwiseAligner, substitution_matrices

from src.download import (
    CURATED_DIR,
    PROJECT_ROOT,
    RAW_DIR,
    UNIPROT_JSON,
    curated_sequence,
    load_seeds,
)

# Curation rule: within a family, no two luciferases may be this similar.
MAX_IDENTITY = 97.0
REQUIRED_COLUMNS = ("label", "organism", "uniprot_url", "reason")
# NCBI records carry no InterPro cross-references.
REQUIRED_UNIPROT_COLUMNS = ("interpro", "taxon_id")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    accession = None
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if accession is not None:
                records.append((accession, "".join(seq)))
            parts = line[1:].split("|")
            if parts[0] in {"sp", "tr"} and len(parts) >= 2:
                accession = parts[1]
            else:
                accession = parts[0].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if accession is not None:
        records.append((accession, "".join(seq)))
    return records


def percent_identity(a: str, b: str) -> float:
    """Global alignment identity, as a percentage of the shorter sequence."""
    aligner = PairwiseAligner(
        mode="global",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )
    aligner.end_gap_score = 0
    alignment = aligner.align(a, b)[0]
    matches = sum(1 for x, y in zip(*alignment) if x == y and x != "-")
    return 100 * matches / min(len(a), len(b))


def verify_local() -> tuple[list[str], list[dict], dict[str, str]]:
    issues: list[str] = []
    seeds = load_seeds()
    expected = [s["accession"] for s in seeds]

    raw_files = {p.name for p in RAW_DIR.glob("*.fasta")}
    if raw_files != {f"{a}.fasta" for a in expected}:
        issues.append(f"raw FASTA set differs from seeds: {sorted(raw_files)}")

    with (CURATED_DIR / "manifest.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [r["accession"] for r in rows] != expected:
        issues.append("manifest accessions/order differ from seeds.csv")
    for seed, row in zip(seeds, rows):
        for col in ("family", "role"):
            if seed[col] != row[col]:
                issues.append(f"{row['accession']} {col} {row[col]!r} != seed {seed[col]!r}")

    seeds_by_acc = {s["accession"]: s for s in seeds}
    seqs: dict[str, str] = {}
    for row in rows:
        acc = row["accession"]
        required = REQUIRED_COLUMNS
        if row["source"] == "uniprot":
            required += REQUIRED_UNIPROT_COLUMNS
        for col in required:
            if not (row.get(col) or "").strip():
                issues.append(f"{acc} empty {col}")
        if row["complete"] != "True":
            issues.append(f"{acc} not marked complete")
        records = parse_fasta((PROJECT_ROOT / row["fasta_path"]).read_text(encoding="utf-8"))
        if len(records) != 1:
            issues.append(f"{acc}: raw FASTA has {len(records)} records")
            continue
        rec_acc, raw_seq = records[0]
        seq = curated_sequence(raw_seq, seeds_by_acc.get(acc, {"accession": acc}))
        seqs[acc] = seq
        if rec_acc != acc:
            issues.append(f"{acc}: raw FASTA accession is {rec_acc}")
        if len(seq) != int(row["length"]):
            issues.append(f"{acc} length manifest {row['length']} vs FASTA {len(seq)}")
        if hashlib.sha256(seq.encode("ascii")).hexdigest() != row["sequence_sha256"]:
            issues.append(f"{acc} sequence hash differs from manifest")
        extra = set(seq) - AMINO_ACIDS
        if extra:
            issues.append(f"{acc} non-standard residues {sorted(extra)}")

    if len(set(seqs.values())) != len(seqs):
        issues.append("identical sequences under different accessions")

    for name, role in (("luciferases.fasta", "luciferase"), ("controls.fasta", "control")):
        combined = parse_fasta((CURATED_DIR / name).read_text(encoding="utf-8"))
        want = [r["accession"] for r in rows if r["role"] == role]
        if [a for a, _ in combined] != want:
            issues.append(f"{name} accessions/order differ from manifest")
        for acc, seq in combined:
            if seqs.get(acc) != seq:
                issues.append(f"{name}: {acc} sequence differs from raw FASTA")
    combined_all = parse_fasta((CURATED_DIR / "all_sequences.fasta").read_text(encoding="utf-8"))
    if [a for a, _ in combined_all] != expected:
        issues.append("all_sequences.fasta accessions/order differ from manifest")

    luciferases = [r for r in rows if r["role"] == "luciferase"]
    for a, b in itertools.combinations(luciferases, 2):
        if a["family"] != b["family"]:
            continue
        ident = percent_identity(seqs[a["accession"]], seqs[b["accession"]])
        if ident >= MAX_IDENTITY:
            issues.append(
                f"near-duplicate luciferases {a['accession']} / {b['accession']}: {ident:.1f}%"
            )

    return issues, rows, seqs


def verify_live(rows: list[dict]) -> list[str]:
    issues: list[str] = []
    session = requests.Session()
    session.headers.update({"User-Agent": "luciferase-course-project/0.2-verify"})
    for i, row in enumerate(rows):
        acc = row["accession"]
        if row["source"] != "uniprot":
            # NCBI accessions are versioned (e.g. .1), so a changed sequence gets a new accession.
            continue
        if i:
            time.sleep(0.35)
        response = session.get(UNIPROT_JSON.format(accession=acc), timeout=60)
        if not response.ok:
            issues.append(f"{acc}: live UniProt returned HTTP {response.status_code}")
            continue
        entry = response.json()
        live_seq = entry.get("sequence", {}).get("value", "")
        if hashlib.sha256(live_seq.encode("ascii")).hexdigest() != row["sequence_sha256"]:
            issues.append(f"{acc}: live UniProt sequence changed since download")
        live_version = str(entry.get("entryAudit", {}).get("sequenceVersion", ""))
        if live_version != row["sequence_version"]:
            issues.append(f"{acc}: sequence version {row['sequence_version']} -> {live_version}")
        organism = entry.get("organism", {}).get("scientificName", "")
        if organism != row["organism"]:
            issues.append(f"{acc}: organism live {organism!r} vs {row['organism']!r}")
        print(f"LIVE {acc} {organism} len={len(live_seq)} type={entry.get('entryType')}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--live", action="store_true", help="also check live UniProt")
    args = parser.parse_args()

    issues, rows, _ = verify_local()
    families: dict[str, dict[str, int]] = {}
    for row in rows:
        families.setdefault(row["family"], {}).setdefault(row["role"], 0)
        families[row["family"]][row["role"]] += 1
    for family, counts in families.items():
        print(f"{family:15s} {counts}")
    if args.live:
        issues += verify_live(rows)

    if issues:
        print("ISSUES")
        for item in issues:
            print("-", item)
        raise SystemExit(1)
    print("ISSUES none")
    print("DATASET VERIFIED")


if __name__ == "__main__":
    main()

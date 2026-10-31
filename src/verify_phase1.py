"""Verify Phase 1 curated dataset against local FASTA and live UniProt."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import requests

from src.download import PROJECT_ROOT, STARTER_RECORDS, fasta_sequence_length

EXPECTED = [r["accession"] for r in STARTER_RECORDS]


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    accession = None
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if accession is not None:
                records.append((accession, "".join(seq)))
            header = line[1:]
            parts = header.split("|")
            if parts and parts[0] in {"sp", "tr"} and len(parts) >= 2:
                accession = parts[1]
            else:
                accession = parts[0].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if accession is not None:
        records.append((accession, "".join(seq)))
    return records


def verify() -> list[str]:
    issues: list[str] = []
    raw = PROJECT_ROOT / "data" / "raw"
    curated = PROJECT_ROOT / "data" / "curated"
    raw_fastas = {p.name for p in raw.glob("*.fasta")}
    expected_files = {f"{a}.fasta" for a in EXPECTED}
    if raw_fastas != expected_files:
        issues.append(f"raw FASTA set mismatch: {sorted(raw_fastas)}")

    with (curated / "manifest.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [r["accession"] for r in rows] != EXPECTED:
        issues.append("manifest accession order mismatch")
    if len(rows) != 11:
        issues.append(f"manifest has {len(rows)} rows, expected 11")

    roles = Counter(r["role"] for r in rows)
    families = Counter(r["family"] for r in rows)
    print("roles", dict(roles))
    print("families", dict(families))
    if roles != Counter({"luciferase": 9, "control": 2}):
        issues.append(f"role counts {roles}")
    expected_families = {
        "firefly": 3,
        "renilla": 2,
        "bacterial_luxa": 2,
        "fungal_luz": 2,
        "control_acs": 1,
        "control_dhla": 1,
    }
    if dict(families) != expected_families:
        issues.append(f"family counts {families}")

    aa = set("ACDEFGHIKLMNPQRSTVWYX")
    for row in rows:
        for col in (
            "protein_name",
            "organism",
            "go_notes",
            "interpro_notes",
            "uniprot_url",
            "reason",
        ):
            if not (row.get(col) or "").strip():
                issues.append(f"{row['accession']} empty {col}")
        if row["complete"] != "True":
            issues.append(f"{row['accession']} not marked complete")
        path = PROJECT_ROOT / row["fasta_path"]
        records = parse_fasta(path.read_text(encoding="utf-8"))
        if len(records) != 1:
            issues.append(f"{path.name} has {len(records)} records")
            continue
        acc, seq = records[0]
        if acc != row["accession"]:
            issues.append(f"{path.name} accession {acc}")
        if len(seq) != int(row["length"]):
            issues.append(
                f"{row['accession']} length manifest {row['length']} vs seq {len(seq)}"
            )
        extra = set(seq) - aa
        if extra:
            issues.append(f"{row['accession']} unexpected residues {extra}")

    luc = parse_fasta((curated / "luciferases.fasta").read_text(encoding="utf-8"))
    ctrl = parse_fasta((curated / "controls.fasta").read_text(encoding="utf-8"))
    print("luciferases n", len(luc), "controls n", len(ctrl))
    if [a for a, _ in luc] != EXPECTED[:9]:
        issues.append(f"luciferases.fasta order {[a for a, _ in luc]}")
    if [a for a, _ in ctrl] != EXPECTED[9:]:
        issues.append(f"controls.fasta order {[a for a, _ in ctrl]}")
    if fasta_sequence_length((curated / "luciferases.fasta").read_text(encoding="utf-8")) != sum(
        int(r["length"]) for r in rows if r["role"] == "luciferase"
    ):
        issues.append("combined luciferases length mismatch")

    session = requests.Session()
    session.headers.update({"User-Agent": "luciferase-course-project/0.1-verify"})
    for row in rows:
        url = f"https://rest.uniprot.org/uniprotkb/{row['accession']}.json"
        entry = session.get(url, timeout=60).json()
        live_len = entry["sequence"]["length"]
        organism = entry["organism"]["scientificName"]
        if live_len != int(row["length"]):
            issues.append(
                f"{row['accession']} live length {live_len} vs {row['length']}"
            )
        if organism != row["organism"]:
            issues.append(
                f"{row['accession']} organism live {organism!r} vs {row['organism']!r}"
            )
        print(
            f"LIVE {row['accession']} {organism} len={live_len} "
            f"type={entry.get('entryType')}"
        )

    return issues


if __name__ == "__main__":
    found = verify()
    if found:
        print("ISSUES")
        for item in found:
            print("-", item)
        raise SystemExit(1)
    print("ISSUES none")
    print("PHASE 1 VERIFIED")

"""Download and curate UniProt luciferase sequences (Phase 1)."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
CURATED_DIR = PROJECT_ROOT / "data" / "curated"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_PAGE = "https://www.uniprot.org/uniprotkb/{accession}"

STARTER_RECORDS = [
    {
        "accession": "P08659",
        "family": "firefly",
        "role": "luciferase",
        "reason": "Classic Photinus pyralis firefly luciferase; reviewed reference sequence.",
    },
    {
        "accession": "P13129",
        "family": "firefly",
        "role": "luciferase",
        "reason": "Nipponoluciola cruciata (formerly Luciola) firefly luciferase.",
    },
    {
        "accession": "Q5UFR2",
        "family": "firefly",
        "role": "luciferase",
        "reason": "Lampyris turkestanicus firefly luciferase; additional beetle representative.",
    },
    {
        "accession": "P27652",
        "family": "renilla",
        "role": "luciferase",
        "reason": "Renilla reniformis luciferase; marine coelenterazine system.",
    },
    {
        "accession": "Q9BLZ3",
        "family": "renilla",
        "role": "luciferase",
        "reason": "Renilla muelleri luciferase; second marine representative.",
    },
    {
        "accession": "P19907",
        "family": "bacterial_luxa",
        "role": "luciferase",
        "reason": "Aliivibrio fischeri LuxA alpha subunit of bacterial luciferase.",
    },
    {
        "accession": "P09140",
        "family": "bacterial_luxa",
        "role": "luciferase",
        "reason": "Photobacterium leiognathi LuxA alpha subunit.",
    },
    {
        "accession": "A0A3G9JYH7",
        "family": "fungal_luz",
        "role": "luciferase",
        "reason": "Neonothopanus nambi fungal Luz luciferase.",
    },
    {
        "accession": "A0A2H3E985",
        "family": "fungal_luz",
        "role": "luciferase",
        "reason": "Armillaria gallica fungal Luz luciferase.",
    },
    {
        "accession": "Q9VCC6",
        "family": "control_acs",
        "role": "control",
        "reason": "Drosophila fatty acyl-CoA synthetase CG6178; non-luminous firefly-family comparison.",
    },
    {
        "accession": "P22643",
        "family": "control_dhla",
        "role": "control",
        "reason": "Xanthobacter haloalkane dehalogenase DhlA; non-luminous Renilla-family comparison.",
    },
]

MANIFEST_FIELDS = [
    "accession",
    "family",
    "role",
    "protein_name",
    "gene",
    "organism",
    "reviewed",
    "length",
    "complete",
    "go_notes",
    "interpro_notes",
    "uniprot_url",
    "fasta_url",
    "reason",
    "fasta_path",
]


def _session() -> requests.Session:
    session = requests.Session()
    session.headers.update({"User-Agent": "luciferase-course-project/0.1 (educational)"})
    return session


def _go_notes(entry: dict, limit: int = 4) -> str:
    terms = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts", []):
            value = (text.get("value") or "").split(".")[0].strip()
            if value:
                terms.append(value)
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        props = {p["key"]: p["value"] for p in xref.get("properties", []) if "key" in p}
        name = props.get("GoTerm", "")
        if name.startswith("F:"):
            terms.append(name[2:])
    uniq = []
    for term in terms:
        if term not in uniq:
            uniq.append(term)
    return "; ".join(uniq[:limit])


def _interpro_notes(entry: dict, limit: int = 4) -> str:
    names = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "InterPro":
            continue
        props = {p["key"]: p["value"] for p in xref.get("properties", []) if "key" in p}
        name = props.get("EntryName")
        acc = xref.get("id")
        if name:
            names.append(f"{acc}:{name}" if acc else name)
    return "; ".join(names[:limit])


def _protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    if rec.get("fullName", {}).get("value"):
        return rec["fullName"]["value"]
    alts = desc.get("submissionNames") or desc.get("alternativeNames") or []
    if alts and alts[0].get("fullName", {}).get("value"):
        return alts[0]["fullName"]["value"]
    return ""


def _gene_name(entry: dict) -> str:
    genes = entry.get("genes") or []
    if not genes:
        return ""
    return (genes[0].get("geneName") or {}).get("value", "")


def fetch_entry(session: requests.Session, accession: str) -> dict:
    url = UNIPROT_JSON.format(accession=accession)
    response = session.get(url, timeout=60)
    response.raise_for_status()
    return response.json()


def fetch_fasta(session: requests.Session, accession: str) -> str:
    url = UNIPROT_FASTA.format(accession=accession)
    response = session.get(url, timeout=60)
    response.raise_for_status()
    text = response.text.strip() + "\n"
    if not text.startswith(">"):
        raise ValueError(f"{accession}: response is not FASTA")
    return text


def fasta_sequence_length(fasta_text: str) -> int:
    lines = fasta_text.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return len(seq.replace(" ", "").replace("*", ""))


def curate_record(session: requests.Session, seed: dict) -> tuple[dict, str]:
    accession = seed["accession"]
    entry = fetch_entry(session, accession)
    fasta = fetch_fasta(session, accession)
    seq_len = fasta_sequence_length(fasta)
    reported = entry.get("sequence", {}).get("length")
    complete = bool(reported) and reported == seq_len and seq_len > 50
    features = entry.get("features") or []
    if any(f.get("type") == "Non-terminal residue" for f in features):
        complete = False

    row = {
        "accession": accession,
        "family": seed["family"],
        "role": seed["role"],
        "protein_name": _protein_name(entry),
        "gene": _gene_name(entry),
        "organism": entry.get("organism", {}).get("scientificName", ""),
        "reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "length": seq_len,
        "complete": complete,
        "go_notes": _go_notes(entry),
        "interpro_notes": _interpro_notes(entry),
        "uniprot_url": UNIPROT_PAGE.format(accession=accession),
        "fasta_url": UNIPROT_FASTA.format(accession=accession),
        "reason": seed["reason"],
        "fasta_path": f"data/raw/{accession}.fasta",
    }
    return row, fasta


def write_combined_fasta(path: Path, items: list[tuple[dict, str]]) -> None:
    chunks = []
    for row, fasta in items:
        header = (
            f">{row['accession']}|{row['family']}|{row['role']}|"
            f"{row['organism'].replace(' ', '_')}|{row['protein_name'].replace(' ', '_')}"
        )
        seq = "".join(
            line.strip() for line in fasta.splitlines() if not line.startswith(">")
        )
        wrapped = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
        chunks.append(f"{header}\n{wrapped}\n")
    path.write_text("".join(chunks), encoding="utf-8")


def run() -> list[dict]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CURATED_DIR.mkdir(parents=True, exist_ok=True)
    session = _session()
    rows: list[dict] = []
    luciferases: list[tuple[dict, str]] = []
    controls: list[tuple[dict, str]] = []

    for i, seed in enumerate(STARTER_RECORDS):
        if i:
            time.sleep(0.35)
        row, fasta = curate_record(session, seed)
        (RAW_DIR / f"{row['accession']}.fasta").write_text(fasta, encoding="utf-8")
        rows.append(row)
        bucket = controls if row["role"] == "control" else luciferases
        bucket.append((row, fasta))
        print(
            f"{row['accession']}\t{row['organism']}\tlen={row['length']}\t"
            f"reviewed={row['reviewed']}\tcomplete={row['complete']}"
        )

    manifest_path = CURATED_DIR / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    write_combined_fasta(CURATED_DIR / "luciferases.fasta", luciferases)
    write_combined_fasta(CURATED_DIR / "controls.fasta", controls)
    (CURATED_DIR / "validation_summary.json").write_text(
        json.dumps(
            {
                "n_records": len(rows),
                "n_luciferases": len(luciferases),
                "n_controls": len(controls),
                "all_complete": all(r["complete"] for r in rows),
                "accessions": [r["accession"] for r in rows],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return rows


if __name__ == "__main__":
    run()

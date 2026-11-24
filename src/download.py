"""Download and curate luciferase and control sequences from UniProt and NCBI.

The sequences to fetch are listed in data/seeds.csv:

- accession, family, role, label (optional display name), reason
- source: "uniprot" (default) or "ncbi" (NCBI Protein, for sequences UniProt lacks)
- allow_partial: "yes" accepts an entry flagged as a fragment/partial when the
  sequence still starts with Met (e.g. CDS clones lacking only a stop codon)
- strip_suffix: residues removed from the C-terminus (e.g. an expression His-tag)
- organism: overrides the record's organism (e.g. for "synthetic construct")

Every record is fetched and checked before any file is written, so a failed run
never leaves a half-updated dataset.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from datetime import date
from pathlib import Path

import requests
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
SEEDS_PATH = DATA_DIR / "seeds.csv"
RAW_DIR = DATA_DIR / "raw"
CURATED_DIR = DATA_DIR / "curated"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_PAGE = "https://www.uniprot.org/uniprotkb/{accession}"
NCBI_EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=protein&id={accession}&rettype={rettype}&retmode=text"
)
NCBI_PAGE = "https://www.ncbi.nlm.nih.gov/protein/{accession}"

FAMILIES = [
    "firefly",
    "renilla_type",
    "bacterial_luxa",
    "fungal_luz",
    "copepod",
    "dinoflagellate",
    "ostracod",
    "click_beetle",
    "photoprotein",
    "decapod_shrimp",
]
ROLES = ["luciferase", "control"]
SOURCES = ["uniprot", "ncbi"]

MANIFEST_FIELDS = [
    "accession",
    "source",
    "family",
    "role",
    "label",
    "uniprot_protein_name",
    "gene",
    "organism",
    "taxon_id",
    "reviewed",
    "evidence_level",
    "length",
    "fragment",
    "partial_accepted",
    "residues_stripped",
    "complete",
    "entry_version",
    "sequence_version",
    "sequence_sha256",
    "function_note",
    "go_molecular_function",
    "interpro",
    "pfam",
    "uniprot_url",
    "fasta_url",
    "reason",
    "fasta_path",
    "retrieved",
]


def load_seeds(path: Path = SEEDS_PATH) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        seeds = list(csv.DictReader(handle))
    seen = set()
    for seed in seeds:
        acc = seed["accession"]
        if acc in seen:
            raise ValueError(f"{acc}: duplicate accession in {path.name}")
        seen.add(acc)
        seed["source"] = seed.get("source") or "uniprot"
        if seed["source"] not in SOURCES:
            raise ValueError(f"{acc}: unknown source {seed['source']!r}")
        if seed["family"] not in FAMILIES:
            raise ValueError(f"{acc}: unknown family {seed['family']!r}")
        if seed["role"] not in ROLES:
            raise ValueError(f"{acc}: unknown role {seed['role']!r}")
    return seeds


def _session() -> requests.Session:
    session = requests.Session()
    session.headers.update({"User-Agent": "luciferase-course-project/0.2 (educational)"})
    return session


def _get(session: requests.Session, url: str, attempts: int = 4) -> requests.Response:
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(url, timeout=60)
            if response.status_code in (429, 500, 502, 503, 504) and attempt < attempts:
                raise requests.HTTPError(f"HTTP {response.status_code}")
            response.raise_for_status()
            return response
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
            if attempt == attempts:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _xref_props(xref: dict) -> dict:
    return {p["key"]: p["value"] for p in xref.get("properties", []) if "key" in p}


def _function_note(entry: dict) -> str:
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts", []):
            value = (text.get("value") or "").split(". ")[0].strip().rstrip(".")
            if value:
                return value
    return ""


def _go_molecular_function(entry: dict, limit: int = 5) -> str:
    terms = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        name = _xref_props(xref).get("GoTerm", "")
        if name.startswith("F:") and name[2:] not in terms:
            terms.append(name[2:])
    return "; ".join(terms[:limit])


def _domain_xrefs(entry: dict, database: str) -> str:
    names = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != database:
            continue
        name = _xref_props(xref).get("EntryName")
        names.append(f"{xref.get('id')}:{name}" if name else xref.get("id", ""))
    return "; ".join(names)


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
    return _get(session, UNIPROT_JSON.format(accession=accession)).json()


def fetch_fasta(session: requests.Session, accession: str) -> str:
    text = _get(session, UNIPROT_FASTA.format(accession=accession)).text.strip() + "\n"
    if not text.startswith(">"):
        raise ValueError(f"{accession}: response is not FASTA")
    return text


def fasta_sequence(fasta_text: str) -> str:
    lines = fasta_text.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq.replace(" ", "").replace("*", "")


def fasta_sequence_length(fasta_text: str) -> int:
    return len(fasta_sequence(fasta_text))


def curated_sequence(raw_seq: str, seed: dict) -> str:
    """Apply the seed's strip_suffix (e.g. a His-tag) to a downloaded sequence."""
    suffix = seed.get("strip_suffix") or ""
    if suffix:
        if not raw_seq.endswith(suffix):
            raise ValueError(f"{seed['accession']}: sequence does not end with {suffix}")
        return raw_seq[: -len(suffix)]
    return raw_seq


def _uniprot_metadata(session: requests.Session, accession: str, raw_seq: str) -> dict:
    entry = fetch_entry(session, accession)
    if entry.get("primaryAccession") != accession:
        raise ValueError(
            f"{accession}: UniProt returned {entry.get('primaryAccession')!r} "
            "(accession merged or demerged; update seeds.csv)"
        )
    if raw_seq != entry.get("sequence", {}).get("value"):
        raise ValueError(f"{accession}: FASTA and JSON sequences differ")
    audit = entry.get("entryAudit", {})
    return {
        "uniprot_protein_name": _protein_name(entry),
        "gene": _gene_name(entry),
        "organism": entry.get("organism", {}).get("scientificName", ""),
        "taxon_id": entry.get("organism", {}).get("taxonId", ""),
        "reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "evidence_level": entry.get("proteinExistence", ""),
        # the flag can also be "Precursor", which is not a fragment
        "fragment": "Fragment" in (entry.get("proteinDescription", {}).get("flag") or ""),
        "nonterminal": any(
            f.get("type") == "Non-terminal residue" for f in entry.get("features") or []
        ),
        "entry_version": audit.get("entryVersion", ""),
        "sequence_version": audit.get("sequenceVersion", ""),
        "function_note": _function_note(entry),
        "go_molecular_function": _go_molecular_function(entry),
        "interpro": _domain_xrefs(entry, "InterPro"),
        "pfam": _domain_xrefs(entry, "Pfam"),
        "uniprot_url": UNIPROT_PAGE.format(accession=accession),
        "fasta_url": UNIPROT_FASTA.format(accession=accession),
    }


def _ncbi_metadata(session: requests.Session, accession: str, raw_seq: str) -> dict:
    gp = _get(session, NCBI_EFETCH.format(accession=accession, rettype="gp")).text
    record = SeqIO.read(io.StringIO(gp), "genbank")
    if str(record.seq) != raw_seq:
        raise ValueError(f"{accession}: GenPept and FASTA sequences differ")
    product, gene, taxon, partial = "", "", "", False
    for feature in record.features:
        q = feature.qualifiers
        if feature.type == "source":
            taxon = next(
                (x.split(":")[1] for x in q.get("db_xref", []) if x.startswith("taxon:")), ""
            )
        elif feature.type == "Protein":
            product = q.get("product", [""])[0]
            partial = "<" in str(feature.location) or ">" in str(feature.location)
        elif feature.type == "CDS":
            gene = q.get("gene", [""])[0]
    version = record.id.rsplit(".", 1)[1] if "." in record.id else ""
    return {
        "uniprot_protein_name": product,
        "gene": gene,
        "organism": record.annotations.get("organism", ""),
        "taxon_id": taxon,
        "reviewed": False,
        "evidence_level": f"NCBI {record.annotations.get('data_file_division', '')} record",
        "fragment": partial or "partial" in record.description.lower(),
        "nonterminal": False,
        "entry_version": version,
        "sequence_version": version,
        "function_note": "",
        "go_molecular_function": "",
        "interpro": "",
        "pfam": "",
        "uniprot_url": NCBI_PAGE.format(accession=accession),
        "fasta_url": NCBI_EFETCH.format(accession=accession, rettype="fasta"),
    }


def fetch_ncbi_fasta(session: requests.Session, accession: str) -> str:
    text = _get(session, NCBI_EFETCH.format(accession=accession, rettype="fasta")).text
    text = text.strip() + "\n"
    if not text.startswith(">"):
        raise ValueError(f"{accession}: NCBI response is not FASTA")
    return text


def curate_record(session: requests.Session, seed: dict) -> tuple[dict, str, str]:
    """Return (manifest row, FASTA as downloaded, curated sequence)."""
    accession = seed["accession"]
    if seed["source"] == "ncbi":
        fasta = fetch_ncbi_fasta(session, accession)
        meta = _ncbi_metadata(session, accession, fasta_sequence(fasta))
    else:
        fasta = fetch_fasta(session, accession)
        meta = _uniprot_metadata(session, accession, fasta_sequence(fasta))
    seq = curated_sequence(fasta_sequence(fasta), seed)

    truncated = meta["fragment"] or meta["nonterminal"]
    # allow_partial also covers full-length mature chains from protein sequencing (no initiator Met)
    no_met = not truncated and not seq.startswith("M")
    partial_accepted = (truncated or no_met) and seed.get("allow_partial") == "yes"
    complete = partial_accepted or (not truncated and seq.startswith("M"))
    row = {
        "accession": accession,
        "source": seed["source"],
        "family": seed["family"],
        "role": seed["role"],
        "label": seed.get("label") or meta["uniprot_protein_name"],
        **{k: v for k, v in meta.items() if k != "nonterminal"},
        "organism": seed.get("organism") or meta["organism"],
        # an overridden organism makes the record's taxon (e.g. synthetic construct) wrong
        "taxon_id": "" if seed.get("organism") else meta["taxon_id"],
        "length": len(seq),
        "partial_accepted": partial_accepted,
        "residues_stripped": seed.get("strip_suffix") or "",
        "complete": complete,
        "sequence_sha256": hashlib.sha256(seq.encode("ascii")).hexdigest(),
        "reason": seed["reason"],
        "fasta_path": f"data/raw/{accession}.fasta",
        "retrieved": date.today().isoformat(),
    }
    return row, fasta, seq


def _header_field(text: str) -> str:
    return text.replace("|", "/").replace(" ", "_")


def write_combined_fasta(path: Path, items: list[tuple[dict, str, str]]) -> None:
    chunks = []
    for row, _, seq in items:
        header = "|".join(
            _header_field(str(row[k]))
            for k in ("accession", "family", "role", "organism", "label")
        )
        wrapped = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
        chunks.append(f">{header}\n{wrapped}\n")
    path.write_text("".join(chunks), encoding="utf-8", newline="\n")


def run() -> list[dict]:
    seeds = load_seeds()
    session = _session()
    items: list[tuple[dict, str, str]] = []

    for i, seed in enumerate(seeds):
        if i:
            time.sleep(0.4)
        items.append(curate_record(session, seed))
        row = items[-1][0]
        print(
            f"{row['accession']}\t{row['family']}/{row['role']}\t{row['organism']}\t"
            f"len={row['length']}\treviewed={row['reviewed']}\tcomplete={row['complete']}"
        )

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CURATED_DIR.mkdir(parents=True, exist_ok=True)
    wanted = {f"{row['accession']}.fasta" for row, _, _ in items}
    for stale in RAW_DIR.glob("*.fasta"):
        if stale.name not in wanted:
            print(f"removing stale {stale.name}")
            stale.unlink()
    for row, fasta, _ in items:
        (RAW_DIR / f"{row['accession']}.fasta").write_text(
            fasta, encoding="utf-8", newline="\n"
        )

    rows = [row for row, _, _ in items]
    with (CURATED_DIR / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    luciferases = [item for item in items if item[0]["role"] == "luciferase"]
    controls = [item for item in items if item[0]["role"] == "control"]
    write_combined_fasta(CURATED_DIR / "luciferases.fasta", luciferases)
    write_combined_fasta(CURATED_DIR / "controls.fasta", controls)
    write_combined_fasta(CURATED_DIR / "all_sequences.fasta", items)

    counts = {
        fam: {
            role: sum(1 for r in rows if r["family"] == fam and r["role"] == role)
            for role in ROLES
        }
        for fam in FAMILIES
    }
    (CURATED_DIR / "validation_summary.json").write_text(
        json.dumps(
            {
                "retrieved": date.today().isoformat(),
                "n_records": len(rows),
                "n_luciferases": len(luciferases),
                "n_controls": len(controls),
                "counts_by_family": counts,
                "all_complete": all(r["complete"] for r in rows),
                "n_reviewed": sum(1 for r in rows if r["reviewed"]),
                "accessions": [r["accession"] for r in rows],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return rows


if __name__ == "__main__":
    run()

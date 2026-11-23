"""Motif search (W2) with multiple pattern matching (W12), plus GO summary (W3).

1. Collect the PROSITE pattern IDs that UniProt lists for our proteins.
2. Download each pattern from PROSITE and convert its syntax to a regex.
3. Scan every protein for every pattern, and record which families each hits.
   A motif shared across independent luciferase families would be evidence of
   sequence-level convergence.
4. Summarise GO molecular-function terms per family from the manifest.

Writes results/prosite_patterns.csv, motif_hits.csv and go_by_family.csv.

Usage:
    python -m src.motifs
"""

from __future__ import annotations

import re
import time

import pandas as pd
import requests

from src.align import load_proteins
from src.download import UNIPROT_JSON
from src.features import RESULTS_DIR

PROSITE_TXT = "https://prosite.expasy.org/{id}.txt"


def prosite_to_regex(pattern: str) -> str:
    """Convert PROSITE pattern syntax, e.g. C-x(2)-[LIVM]-{P}, to a Python regex."""
    regex = pattern.rstrip(".").replace("<", "^").replace(">", "$")
    regex = regex.replace("{", "[^").replace("}", "]").replace("(", "{").replace(")", "}")
    return regex.replace("-", "").replace("x", ".")


def prosite_ids(proteins: pd.DataFrame) -> set[str]:
    ids, session = set(), requests.Session()
    for acc in proteins.loc[proteins["source"] == "uniprot", "accession"]:
        entry = session.get(UNIPROT_JSON.format(accession=acc), timeout=60).json()
        ids |= {x["id"] for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "PROSITE"}
        time.sleep(0.2)
    return ids


def fetch_patterns(ids: set[str]) -> pd.DataFrame:
    rows = []
    for pid in sorted(ids):
        text = requests.get(PROSITE_TXT.format(id=pid), timeout=60).text
        fields = {line[:2]: line[5:].strip() for line in text.splitlines() if line[:2] in ("ID", "DE")}
        pattern = "".join(line[5:].strip() for line in text.splitlines() if line.startswith("PA"))
        if pattern:  # profiles (no PA line) cannot be expressed as a regex
            rows.append({"prosite_id": pid, "name": fields.get("ID", "").split(";")[0],
                         "description": fields.get("DE", ""), "pattern": pattern,
                         "regex": prosite_to_regex(pattern)})
        time.sleep(0.2)
    return pd.DataFrame(rows)


def scan(proteins: pd.DataFrame, patterns: pd.DataFrame) -> pd.DataFrame:
    compiled = [(p.prosite_id, p.name, re.compile(p.regex)) for p in patterns.itertuples()]
    rows = []
    for prot in proteins.itertuples():
        for pid, name, regex in compiled:
            for m in regex.finditer(prot.sequence):
                rows.append({"prosite_id": pid, "name": name, "accession": prot.accession,
                             "family": prot.family, "role": prot.role, "start": m.start() + 1,
                             "match": m.group()})
    return pd.DataFrame(rows)


def go_by_family(proteins: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (family, role), group in proteins.groupby(["family", "role"], observed=True, sort=False):
        terms = sorted({t.strip() for v in group["go_molecular_function"] for t in v.split(";") if t.strip()})
        rows.append({"family": family, "role": role, "n": len(group), "go_molecular_function": "; ".join(terms)})
    return pd.DataFrame(rows)


def run() -> pd.DataFrame:
    proteins = load_proteins()
    patterns = fetch_patterns(prosite_ids(proteins))
    patterns.to_csv(RESULTS_DIR / "prosite_patterns.csv", index=False, lineterminator="\n")
    hits = scan(proteins, patterns)
    hits.to_csv(RESULTS_DIR / "motif_hits.csv", index=False, lineterminator="\n")
    go_by_family(proteins).to_csv(RESULTS_DIR / "go_by_family.csv", index=False, lineterminator="\n")
    luc = hits[hits["role"] == "luciferase"]
    summary = luc.groupby(["prosite_id", "name"])["family"].agg(lambda f: ", ".join(sorted(set(map(str, f)))))
    print(f"{len(patterns)} PROSITE patterns; {len(hits)} matches")
    print(summary.to_string())
    return hits


if __name__ == "__main__":
    run()

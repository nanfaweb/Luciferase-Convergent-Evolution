"""All-vs-all protein BLAST over the curated dataset with local BLAST+.

BLAST+ cannot handle spaces in file paths (it expands relative paths too), and
this project's folder has one. So BLAST runs in a throwaway system temp folder.

Writes results/blast_all_vs_all.tsv and results/blast_summary.csv.

Usage:
    python -m src.blast
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from src.download import CURATED_DIR
from src.features import RESULTS_DIR

EVALUE = 1e-3
FIELDS = "qseqid sseqid pident length evalue bitscore"
FALLBACK_BIN = Path.home() / "tools" / "ncbi-blast-2.17.0+" / "bin"


def blast_tool(name: str) -> str:
    found = shutil.which(name) or shutil.which(str(FALLBACK_BIN / name))
    if not found:
        raise FileNotFoundError(f"{name} not found; install BLAST+ (see CLAUDE.md)")
    return found


def _run(args: list[str], cwd: str) -> str:
    result = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{Path(args[0]).name} failed: {result.stderr.strip()}")
    return result.stdout


def run_blast() -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="luciferase_blast_") as work:
        if " " in work:
            raise RuntimeError(f"temp folder {work!r} has a space; BLAST+ cannot use it")
        shutil.copy(CURATED_DIR / "all_sequences.fasta", Path(work) / "all.fa")
        _run([blast_tool("makeblastdb"), "-in", "all.fa", "-dbtype", "prot", "-out", "allseq"], work)
        stdout = _run(
            [blast_tool("blastp"), "-query", "all.fa", "-db", "allseq", "-evalue", str(EVALUE),
             "-max_target_seqs", "100", "-outfmt", f"6 {FIELDS}"],
            work,
        )
    rows = [line.split("\t") for line in stdout.splitlines() if line]
    hits = pd.DataFrame(rows, columns=FIELDS.split())
    for col in ("query", "subject"):
        parts = hits[f"{col[0]}seqid"].str.split("|", expand=True)
        hits[f"{col}_accession"], hits[f"{col}_family"], hits[f"{col}_role"] = parts[0], parts[1], parts[2]
    hits = hits.drop(columns=["qseqid", "sseqid"])
    hits[["pident", "evalue", "bitscore"]] = hits[["pident", "evalue", "bitscore"]].astype(float)
    # keep the best HSP per query-subject pair, and drop self hits
    hits = (hits.sort_values("bitscore", ascending=False)
            .drop_duplicates(["query_accession", "subject_accession"]))
    return hits[hits["query_accession"] != hits["subject_accession"]].reset_index(drop=True)


def summarise(hits: pd.DataFrame) -> pd.DataFrame:
    hits = hits.assign(same_family=hits["query_family"] == hits["subject_family"])
    manifest = pd.read_csv(CURATED_DIR / "manifest.csv", keep_default_na=False)
    rows = []
    for q in manifest.itertuples():
        own = hits[hits["query_accession"] == q.accession]
        cross = own[~own["same_family"]]
        rows.append({
            "accession": q.accession, "family": q.family, "role": q.role,
            "same_family_hits": int(own["same_family"].sum()),
            "cross_family_hits": len(cross),
            "cross_family_subjects": " ".join(cross["subject_accession"]),
            "best_cross_family_evalue": cross["evalue"].min() if len(cross) else None,
        })
    return pd.DataFrame(rows)


def run() -> pd.DataFrame:
    hits = run_blast()
    hits.to_csv(RESULTS_DIR / "blast_all_vs_all.tsv", sep="\t", index=False, lineterminator="\n")
    summary = summarise(hits)
    summary.to_csv(RESULTS_DIR / "blast_summary.csv", index=False, lineterminator="\n")
    print(f"{len(hits)} hits at E < {EVALUE}; "
          f"{int((summary['cross_family_hits'] > 0).sum())} of {len(summary)} queries hit another family")
    print(summary.groupby(["family", "role"], sort=False)[["same_family_hits", "cross_family_hits"]].sum())
    return summary


if __name__ == "__main__":
    run()

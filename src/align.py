"""Pairwise sequence similarity for every pair of curated proteins.

- Global (Needleman-Wunsch) alignment: percent identity over the shorter sequence.
- Local (Smith-Waterman) alignment: raw score, plus a z-score against shuffled
  copies of the second sequence, which shows whether a score beats chance.

Writes results/pairwise.csv, identity_matrix.csv, similarity_summary.csv,
example_alignments.txt and figures/identity_heatmap.png.

Usage:
    python -m src.align
"""

from __future__ import annotations

import itertools
import random

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

from src.download import CURATED_DIR, FAMILIES
from src.features import RESULTS_DIR
from src.verify_data import parse_fasta

N_SHUFFLES = 20
EXAMPLE_PAIRS = [("P08659", "P13129"), ("P08659", "P19907")]


def make_aligner(mode: str) -> PairwiseAligner:
    aligner = PairwiseAligner(
        mode=mode,
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )
    if mode == "global":
        aligner.end_gap_score = 0
    return aligner


GLOBAL = make_aligner("global")
LOCAL = make_aligner("local")


def global_identity(a: str, b: str) -> float:
    alignment = GLOBAL.align(a, b)[0]
    matches = sum(1 for x, y in zip(*alignment) if x == y and x != "-")
    return 100 * matches / min(len(a), len(b))


def local_zscore(a: str, b: str, rng: random.Random) -> tuple[float, float]:
    score = LOCAL.score(a, b)
    shuffled = []
    for _ in range(N_SHUFFLES):
        letters = list(b)
        rng.shuffle(letters)
        shuffled.append(LOCAL.score(a, "".join(letters)))
    sd = np.std(shuffled)
    return score, (score - np.mean(shuffled)) / sd if sd else 0.0


def load_proteins() -> pd.DataFrame:
    manifest = pd.read_csv(CURATED_DIR / "manifest.csv", keep_default_na=False)
    seqs = dict(parse_fasta((CURATED_DIR / "all_sequences.fasta").read_text(encoding="utf-8")))
    manifest["sequence"] = manifest["accession"].map(seqs)
    manifest["family"] = pd.Categorical(manifest["family"], FAMILIES, ordered=True)
    # luciferases first within each family, then that family's controls
    return manifest.sort_values(
        ["family", "role"], ascending=[True, False], kind="stable"
    ).reset_index(drop=True)


def pair_category(a: pd.Series, b: pd.Series) -> str:
    if a["family"] != b["family"]:
        return "cross-family" if a["role"] == b["role"] == "luciferase" else "cross-family (control)"
    if a["role"] == b["role"] == "luciferase":
        return "same family"
    if a["role"] != b["role"]:
        return "luciferase vs own control"
    return "control vs control"


def pairwise_table(proteins: pd.DataFrame) -> pd.DataFrame:
    rng = random.Random(1)
    rows = []
    for i, j in itertools.combinations(range(len(proteins)), 2):
        a, b = proteins.iloc[i], proteins.iloc[j]
        sw_score, sw_z = local_zscore(a["sequence"], b["sequence"], rng)
        rows.append(
            {
                "accession_a": a["accession"],
                "accession_b": b["accession"],
                "family_a": a["family"],
                "family_b": b["family"],
                "category": pair_category(a, b),
                "global_identity_pct": global_identity(a["sequence"], b["sequence"]),
                "local_score": sw_score,
                "local_zscore": sw_z,
            }
        )
    return pd.DataFrame(rows)


def identity_matrix(proteins: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    acc = list(proteins["accession"])
    matrix = pd.DataFrame(100.0, index=acc, columns=acc)
    for row in pairs.itertuples():
        matrix.loc[row.accession_a, row.accession_b] = row.global_identity_pct
        matrix.loc[row.accession_b, row.accession_a] = row.global_identity_pct
    return matrix


def plot_heatmap(proteins: pd.DataFrame, matrix: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11, 10))
    image = ax.imshow(matrix.values, cmap="viridis", vmin=0, vmax=100)
    labels = [
        f"{r.accession} {r.family}{' (ctrl)' if r.role == 'control' else ''}"
        for r in proteins.itertuples()
    ]
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=5)
    ax.set_yticks(range(len(labels)), labels, fontsize=5)
    # white lines between families
    bounds = np.flatnonzero(proteins["family"].astype(str).ne(proteins["family"].astype(str).shift()))[1:]
    for b in bounds:
        ax.axhline(b - 0.5, color="white", lw=0.8)
        ax.axvline(b - 0.5, color="white", lw=0.8)
    fig.colorbar(image, ax=ax, shrink=0.7, label="Global identity (%)")
    ax.set_title("Pairwise identity (Needleman-Wunsch, BLOSUM62)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "identity_heatmap.png", dpi=200)
    plt.close(fig)


def write_examples(proteins: pd.DataFrame) -> None:
    seqs = dict(zip(proteins["accession"], proteins["sequence"]))
    fam = dict(zip(proteins["accession"], proteins["family"].astype(str)))
    out = []
    for a, b in EXAMPLE_PAIRS:
        for name, aligner in (("Global (Needleman-Wunsch)", GLOBAL), ("Local (Smith-Waterman)", LOCAL)):
            alignment = aligner.align(seqs[a], seqs[b])[0]
            out.append(
                f"=== {name}: {a} ({fam[a]}) vs {b} ({fam[b]}) | score {alignment.score:.1f} ===\n"
                f"{alignment}\n"
            )
    (RESULTS_DIR / "example_alignments.txt").write_text("\n".join(out), encoding="utf-8", newline="\n")


DP_EXAMPLE = ("HEAGAWGHEE", "PAWHEAE")  # the classic textbook pair (Durbin et al.)
DP_GAP = -8


def needleman_wunsch_dp(a: str, b: str, gap: int = DP_GAP) -> tuple[np.ndarray, str, str]:
    """Textbook dynamic programming: fill the score table, then trace back."""
    blosum = substitution_matrices.load("BLOSUM62")
    table = np.zeros((len(a) + 1, len(b) + 1), dtype=int)
    table[:, 0] = gap * np.arange(len(a) + 1)
    table[0, :] = gap * np.arange(len(b) + 1)
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            table[i, j] = max(table[i - 1, j - 1] + blosum[a[i - 1], b[j - 1]],
                              table[i - 1, j] + gap, table[i, j - 1] + gap)
    top, bottom, i, j = [], [], len(a), len(b)
    while i or j:
        if i and j and table[i, j] == table[i - 1, j - 1] + blosum[a[i - 1], b[j - 1]]:
            top.append(a[i - 1]); bottom.append(b[j - 1]); i, j = i - 1, j - 1
        elif i and table[i, j] == table[i - 1, j] + gap:
            top.append(a[i - 1]); bottom.append("-"); i -= 1
        else:
            top.append("-"); bottom.append(b[j - 1]); j -= 1
    return table, "".join(reversed(top)), "".join(reversed(bottom))


def write_dp_example() -> None:
    a, b = DP_EXAMPLE
    table, top, bottom = needleman_wunsch_dp(a, b)
    check = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                            open_gap_score=DP_GAP, extend_gap_score=DP_GAP).score(a, b)
    header = "      -" + "".join(f"{c:>5}" for c in b)
    lines = [f"Needleman-Wunsch DP table: {a} vs {b} (BLOSUM62, linear gap {DP_GAP})", "", header]
    for i, row in enumerate(table):
        lines.append(f"{('-' + a)[i]:>2}" + "".join(f"{v:>5}" for v in row))
    lines += ["", f"Optimal alignment (score {table[-1, -1]}):", top, bottom, "",
              f"Biopython score for the same settings: {check:.0f} "
              f"({'matches' if check == table[-1, -1] else 'DOES NOT match'} the hand-written DP)"]
    (RESULTS_DIR / "dp_example.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print("\n".join(lines[-4:]))


def run() -> pd.DataFrame:
    (RESULTS_DIR / "figures").mkdir(parents=True, exist_ok=True)
    write_dp_example()
    proteins = load_proteins()
    pairs = pairwise_table(proteins)
    pairs.to_csv(RESULTS_DIR / "pairwise.csv", index=False, float_format="%.3f", lineterminator="\n")
    matrix = identity_matrix(proteins, pairs)
    matrix.to_csv(RESULTS_DIR / "identity_matrix.csv", float_format="%.1f", lineterminator="\n")
    summary = pairs.groupby("category")[["global_identity_pct", "local_zscore"]].agg(
        ["count", "median", "min", "max"]
    )
    summary.columns = [f"{c}_{s}" for c, s in summary.columns]
    summary.to_csv(RESULTS_DIR / "similarity_summary.csv", float_format="%.2f", lineterminator="\n")
    plot_heatmap(proteins, matrix)
    write_examples(proteins)
    print(summary.to_string())
    return pairs


if __name__ == "__main__":
    run()

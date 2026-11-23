"""Dot plots: a picture of where two sequences match.

A dot is drawn at (i, j) when a window of 10 residues starting at position i of
sequence A and position j of sequence B scores at least THRESHOLD with BLOSUM62.
Related sequences show a diagonal line; unrelated ones show scattered noise.

Writes results/figures/dotplots.png.

Usage:
    python -m src.dotplot
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align import substitution_matrices

from src.align import load_proteins
from src.features import RESULTS_DIR

WINDOW = 10
THRESHOLD = 23  # about 2.3 BLOSUM62 points per residue
PANELS = [
    ("P08659", "P13129", "Same family: two firefly luciferases"),
    ("P08659", "A9ZPM4", "Firefly luciferase vs non-luminous control"),
    ("P08659", "P19907", "Different families: firefly vs bacterial LuxA"),
]


def window_scores(a: str, b: str) -> np.ndarray:
    blosum = substitution_matrices.load("BLOSUM62")
    alphabet = blosum.alphabet
    ia = np.array([alphabet.index(c) for c in a])
    ib = np.array([alphabet.index(c) for c in b])
    scores = np.asarray(blosum)[np.ix_(ia, ib)]
    n, m = len(a) - WINDOW + 1, len(b) - WINDOW + 1
    total = np.zeros((n, m))
    for k in range(WINDOW):  # sum along each diagonal window
        total += scores[k : k + n, k : k + m]
    return total


def run() -> None:
    proteins = load_proteins().set_index("accession")
    fig, axes = plt.subplots(1, len(PANELS), figsize=(15, 5))
    for ax, (a, b, title) in zip(axes, PANELS):
        hits = np.argwhere(window_scores(proteins.at[a, "sequence"], proteins.at[b, "sequence"]) >= THRESHOLD)
        ax.scatter(hits[:, 1], hits[:, 0], s=0.4, color="black")
        ax.set_xlim(0, len(proteins.at[b, "sequence"]))
        ax.set_ylim(len(proteins.at[a, "sequence"]), 0)
        ax.set_xlabel(f"{b} ({proteins.at[b, 'organism']})")
        ax.set_ylabel(f"{a} ({proteins.at[a, 'organism']})")
        ax.set_title(title, fontsize=10)
        print(f"{a} vs {b}: {len(hits)} dots")
    fig.suptitle(f"Dot plots (window {WINDOW}, BLOSUM62 score >= {THRESHOLD})")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "dotplots.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    run()

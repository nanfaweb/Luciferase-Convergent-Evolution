"""Within-family neighbour-joining trees from pairwise global identity.

Distance is the Poisson-corrected protein distance d = -ln(1 - p), where p is
the fraction of differing residues from results/identity_matrix.csv. Each tree
is rooted on the family's control most distant from its luciferases, or at the
midpoint when the family has no control. Families are never mixed in one tree.

Writes results/trees/<family>.nwk and results/figures/trees.png.

Usage:
    python -m src.phylo   (run src.align first)
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from src.align import load_proteins
from src.download import FAMILIES
from src.features import RESULTS_DIR


def poisson_distance(identity_pct: float) -> float:
    p = min(1 - identity_pct / 100, 0.99)
    return -math.log(1 - p)


def family_tree(members: pd.DataFrame, identity: pd.DataFrame):
    names = list(members["accession"])
    matrix = [[poisson_distance(identity.loc[a, b]) for b in names[: i + 1]] for i, a in enumerate(names)]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, matrix))
    controls = members.loc[members["role"] == "control", "accession"].tolist()
    luciferases = members.loc[members["role"] == "luciferase", "accession"].tolist()
    if controls:
        # the control least identical to the family's luciferases
        outgroup = min(controls, key=lambda c: identity.loc[c, luciferases].sum())
        tree.root_with_outgroup(outgroup)
    else:
        tree.root_at_midpoint()
    labels = dict(zip(members["accession"], members["organism"]))
    roles = dict(zip(members["accession"], members["role"]))
    for clade in tree.get_terminals():
        acc = clade.name
        suffix = " [control]" if roles[acc] == "control" else ""
        clade.name = f"{acc} {' '.join(labels[acc].split()[:2])}{suffix}"
    return tree


def run() -> None:
    proteins = load_proteins()
    identity = pd.read_csv(RESULTS_DIR / "identity_matrix.csv", index_col=0)
    (RESULTS_DIR / "trees").mkdir(exist_ok=True)
    rows = math.ceil(len(FAMILIES) / 2)
    fig, axes = plt.subplots(rows, 2, figsize=(14, 5 * rows))
    for ax, family in zip(axes.flat, FAMILIES):
        tree = family_tree(proteins[proteins["family"] == family], identity)
        Phylo.write(tree, RESULTS_DIR / "trees" / f"{family}.nwk", "newick")
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda c: c.name if c.is_terminal() else None)
        ax.set_title(family)
        ax.set_ylabel("")
        ax.set_xlabel("Poisson distance (substitutions per site)")
        print(f"{family}: {tree.count_terminals()} sequences")
    for ax in axes.flat[len(FAMILIES):]:
        ax.axis("off")
    fig.suptitle("Within-family neighbour-joining trees (no cross-family tree by design)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(RESULTS_DIR / "figures" / "trees.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()

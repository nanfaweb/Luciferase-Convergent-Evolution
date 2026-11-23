"""Clustering and statistics on the protein properties.

1. Clustering: standardized composition-type features (GRAVY, pI, aromaticity,
   20 amino-acid percentages; not length or weight, which mostly reflect
   protein architecture). PCA plot, and Ward clustering compared with the
   family labels by the adjusted Rand index (1 = groups exactly by family,
   0 = no better than chance).
2. Kruskal-Wallis test across all luciferase families for each main feature.
3. Luciferase vs own control: direction of the median difference in each
   family that has controls, with a two-sided sign test.

Writes results/stats.csv, control_comparison.csv,
figures/pca_features.png and figures/feature_boxplots.png.

Usage:
    python -m src.analysis   (run src.features first)
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import binomtest, kruskal
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

from src.download import FAMILIES
from src.features import RESULTS_DIR

MAIN_FEATURES = ["gravy", "isoelectric_point", "aromaticity", "aa_C_pct"]


def composition_columns(table: pd.DataFrame) -> list[str]:
    return ["gravy", "isoelectric_point", "aromaticity"] + [c for c in table if c.startswith("aa_")]


def clustering(table: pd.DataFrame) -> dict:
    luc = table[table["role"] == "luciferase"].reset_index(drop=True)
    x = StandardScaler().fit_transform(luc[composition_columns(luc)])
    clusters = fcluster(linkage(x, "ward"), t=len(FAMILIES), criterion="maxclust")
    ari = adjusted_rand_score(luc["family"], clusters)

    pca = PCA(n_components=2).fit(x)
    coords = pca.transform(x)
    fig, ax = plt.subplots(figsize=(8, 6))
    for family in FAMILIES:
        mask = luc["family"] == family
        ax.scatter(coords[mask, 0], coords[mask, 1], label=family, s=30)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.0%} of variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.0%} of variance)")
    ax.set_title(f"Luciferases by composition features (Ward clusters vs family: ARI = {ari:.2f})")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "pca_features.png", dpi=200)
    plt.close(fig)
    return {"test": "Ward clustering vs family (adjusted Rand index)", "feature": "all composition",
            "statistic": ari, "p_value": np.nan, "n": len(luc)}


def kruskal_tests(table: pd.DataFrame) -> list[dict]:
    luc = table[table["role"] == "luciferase"]
    rows = []
    for feature in MAIN_FEATURES:
        groups = [luc.loc[luc["family"] == f, feature] for f in FAMILIES]
        stat, p = kruskal(*groups)
        rows.append({"test": f"Kruskal-Wallis across {len(FAMILIES)} luciferase families", "feature": feature,
                     "statistic": stat, "p_value": p, "n": len(luc)})
    return rows


def control_comparison(table: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    rows = []
    for family in FAMILIES:
        fam = table[table["family"] == family]
        if not (fam["role"] == "control").any():
            continue
        for feature in MAIN_FEATURES:
            luc = fam.loc[fam["role"] == "luciferase", feature].median()
            ctrl = fam.loc[fam["role"] == "control", feature].median()
            rows.append({"family": family, "feature": feature, "luciferase_median": luc,
                         "control_median": ctrl, "difference": luc - ctrl,
                         "direction": "higher" if luc > ctrl else "lower"})
    comparison = pd.DataFrame(rows)
    tests = []
    for feature, group in comparison.groupby("feature", sort=False):
        higher = int((group["direction"] == "higher").sum())
        tests.append({"test": f"Sign test: luciferase vs own control ({len(group)} families)", "feature": feature,
                      "statistic": higher, "p_value": binomtest(higher, len(group)).pvalue,
                      "n": len(group)})
    return comparison, tests


def boxplots(table: pd.DataFrame) -> None:
    labels = {"gravy": "GRAVY (hydrophobicity)", "isoelectric_point": "Isoelectric point (pI)",
              "aromaticity": "Aromaticity", "aa_C_pct": "Cysteine (%)"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, feature in zip(axes.flat, MAIN_FEATURES):
        luc = [table.loc[(table["family"] == f) & (table["role"] == "luciferase"), feature] for f in FAMILIES]
        ax.boxplot(luc, tick_labels=FAMILIES)
        for i, f in enumerate(FAMILIES, start=1):
            ctrl = table.loc[(table["family"] == f) & (table["role"] == "control"), feature]
            ax.scatter([i] * len(ctrl), ctrl, marker="x", color="red", zorder=3,
                       label="control" if i == 1 else None)
        ax.set_title(labels[feature])
        ax.tick_params(axis="x", rotation=30, labelsize=8)
    axes.flat[0].legend()
    fig.suptitle("Luciferase properties by family (boxes) with non-luminous controls (red x)")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "feature_boxplots.png", dpi=200)
    plt.close(fig)


def run() -> pd.DataFrame:
    table = pd.read_csv(RESULTS_DIR / "features.csv")
    comparison, sign_tests = control_comparison(table)
    stats = pd.DataFrame([clustering(table)] + kruskal_tests(table) + sign_tests)
    stats.to_csv(RESULTS_DIR / "stats.csv", index=False, float_format="%.4g", lineterminator="\n")
    comparison.to_csv(RESULTS_DIR / "control_comparison.csv", index=False, float_format="%.4f",
                      lineterminator="\n")
    boxplots(table)
    print(stats.to_string(index=False))
    print(comparison.pivot(index="feature", columns="family", values="direction").to_string())
    return stats


if __name__ == "__main__":
    run()

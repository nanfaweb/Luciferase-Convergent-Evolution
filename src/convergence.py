"""Novelty phase: a pre-registered, family-aware test of physicochemical convergence.

The earlier observation was that luciferases have lower GRAVY and higher aromaticity than
their own non-luminous relatives in all 5 original families with controls. The prediction
was written into PLAN.md (git commit d2978e2) *before* two new, independent luminous
families were downloaded:

- photoprotein: Ca2+-regulated photoproteins vs EF-hand SCP/calmodulin;
- decapod_shrimp: Oplophorus luciferase vs FABPs.

This module tests:
1. Replication: do both new families show both predicted directions?
2. A one-sided sign test over all families with controls (the unit is the family).
3. A within-family label-permutation test. Labels are shuffled only inside each family, so
   family ancestry cannot create the signal. Uses every protein. 10,000 permutations,
   seed 1, Holm correction over the 2 predicted properties.
4. Robustness: firefly + click beetle merged into one ACS origin; photoproteins excluded;
   only the original 5 families; leave-one-family-out.
5. Calibration: 500 random features should give p < 0.05 about 5% of the time.
6. Exploratory (not pre-registered): substrate convergence. Are the 4 unrelated families
   that use coelenterazine closer in composition than random sets of 4 of the 10 families?
   Exact enumeration of all 210 label sets.

Writes results/convergence_tests.csv, convergence_by_family.csv, substrate_convergence.csv,
figures/h2_replication.png and figures/substrate_distance.png.

Usage:
    python -m src.convergence   (needs src.features output)
"""

from __future__ import annotations

import itertools

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from src.download import FAMILIES
from src.features import RESULTS_DIR

PREDICTED = {"gravy": -1, "aromaticity": +1}      # sign of (luciferase - control), fixed in advance
EXPLORATORY = ["isoelectric_point", "aa_C_pct"]
ORIGINAL = ["firefly", "renilla_type", "bacterial_luxa", "fungal_luz", "click_beetle"]
NEW = ["photoprotein", "decapod_shrimp"]
N_PERM = 10_000
LUCIFERIN = {
    "firefly": "D-luciferin", "click_beetle": "D-luciferin",
    "renilla_type": "coelenterazine", "copepod": "coelenterazine",
    "photoprotein": "coelenterazine", "decapod_shrimp": "coelenterazine",
    "bacterial_luxa": "reduced flavin + aldehyde", "fungal_luz": "3-hydroxyhispidin",
    "dinoflagellate": "tetrapyrrole", "ostracod": "Cypridina luciferin",
}
NAMES = {"firefly": "Firefly", "renilla_type": "Renilla-type", "bacterial_luxa": "Bacterial LuxA",
         "fungal_luz": "Fungal Luz", "copepod": "Copepod", "dinoflagellate": "Dinoflagellate",
         "ostracod": "Ostracod", "click_beetle": "Click beetle", "photoprotein": "Photoprotein",
         "decapod_shrimp": "Oplophorus shrimp"}


def load() -> pd.DataFrame:
    table = pd.read_csv(RESULTS_DIR / "features.csv", keep_default_na=False)
    for col in [c for c in table if c.startswith("aa_")] + ["gravy", "aromaticity", "isoelectric_point"]:
        table[col] = table[col].astype(float)
    return table


def with_controls(table: pd.DataFrame) -> list[str]:
    roles = table.groupby("family")["role"].agg(set)
    return [f for f in FAMILIES if f in roles and roles[f] == {"luciferase", "control"}]


def family_differences(table: pd.DataFrame, families: list[str]) -> pd.DataFrame:
    """Median(luciferase) - median(control), raw and in units of the feature's overall SD."""
    rows = []
    for feature in list(PREDICTED) + EXPLORATORY:
        sd = table[feature].std()
        for fam in families:
            g = table[table["family"] == fam]
            diff = g.loc[g.role == "luciferase", feature].median() - g.loc[g.role == "control", feature].median()
            rows.append({"family": fam, "feature": feature, "difference": diff, "difference_sd": diff / sd,
                         "predicted_direction": {-1: "lower", 1: "higher"}.get(PREDICTED.get(feature), ""),
                         "observed_direction": "higher" if diff > 0 else "lower",
                         "set": "new (pre-registered)" if fam in NEW else "original"})
    return pd.DataFrame(rows)


def permutation_p(table: pd.DataFrame, feature: str, families: list[str], direction: int | None,
                  n_perm: int = N_PERM, seed: int = 1, groups: dict | None = None) -> tuple[float, float]:
    """Within-family label permutation. Statistic: sum over families of
    (mean luciferase z - mean control z), times the predicted direction (one-sided).
    direction=None gives a two-sided test. `groups` can merge families into one stratum."""
    rng = np.random.default_rng(seed)
    z = (table[feature] - table[feature].mean()) / table[feature].std()
    strata = table["family"].map(groups or {}).fillna(table["family"])
    keep = table["family"].isin(families)
    observed, null = 0.0, np.zeros(n_perm)
    for _, idx in table[keep].groupby(strata[keep]).groups.items():
        values = z[idx].to_numpy()
        is_luc = (table.loc[idx, "role"] == "luciferase").to_numpy()
        n_l, n_c, total = is_luc.sum(), (~is_luc).sum(), values.sum()
        if n_l == 0 or n_c == 0:
            continue
        s = values[is_luc].sum()
        observed += s / n_l - (total - s) / n_c
        perm_s = rng.permuted(np.tile(values, (n_perm, 1)), axis=1)[:, :n_l].sum(axis=1)
        null += perm_s / n_l - (total - perm_s) / n_c
    if direction is None:
        p = (np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1)
    else:
        p = (np.sum(direction * null >= direction * observed) + 1) / (n_perm + 1)
    return observed, p


def sign_test(diffs: pd.DataFrame, feature: str, families: list[str]) -> tuple[int, int, float]:
    d = diffs[(diffs.feature == feature) & diffs.family.isin(families)]["difference"]
    k = int((np.sign(d) == PREDICTED[feature]).sum())
    return k, len(d), binomtest(k, len(d), 0.5, alternative="greater").pvalue


def holm(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    adjusted, running = np.empty(len(pvalues)), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[i]))
        adjusted[i] = running
    return list(adjusted)


def run_tests(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    families = with_controls(table)
    diffs = family_differences(table, families)
    rows = []

    def add(analysis, fams, groups=None, label_families=None):
        ps = []
        for feature, direction in PREDICTED.items():
            k, n, p_sign = sign_test(diffs, feature, label_families or fams)
            stat, p_perm = permutation_p(table, feature, fams, direction, groups=groups)
            ps.append(p_perm)
            rows.append({"analysis": analysis, "feature": feature, "families": n,
                         "families_in_predicted_direction": k, "sign_test_p_one_sided": p_sign,
                         "permutation_statistic": stat, "permutation_p_one_sided": p_perm})
        for row, adj in zip(rows[-len(PREDICTED):], holm(ps)):
            row["permutation_p_holm"] = adj

    # 1. pre-registered replication (new families only), 2-3. main test on all families
    add("Replication: new families only (pre-registered)", NEW)
    add("Main: all families with controls", families)
    # 4. robustness
    add("Robustness: firefly + click beetle as one ACS origin", families,
        groups={"click_beetle": "firefly"}, label_families=[f for f in families if f != "click_beetle"])
    add("Robustness: photoproteins excluded", [f for f in families if f != "photoprotein"])
    add("Exploratory: original 5 families only (hint, before replication)", ORIGINAL)
    for left_out in families:
        add(f"Leave-one-out: without {left_out}", [f for f in families if f != left_out])
    # non-predicted features, two-sided
    for feature in EXPLORATORY:
        stat, p = permutation_p(table, feature, families, None)
        rows.append({"analysis": "Not predicted (two-sided): all families", "feature": feature,
                     "families": len(families), "permutation_statistic": stat, "permutation_p_two_sided": p})
    # 5. calibration with random features
    rng = np.random.default_rng(7)
    fake = table.copy()
    hits = 0
    for i in range(500):
        fake["random"] = rng.normal(size=len(fake))
        _, p = permutation_p(fake, "random", families, 1, n_perm=999, seed=i)
        hits += p < 0.05
    rows.append({"analysis": "Calibration: 500 random features (expect about 0.05)", "feature": "random",
                 "families": len(families), "permutation_p_one_sided": hits / 500})
    return pd.DataFrame(rows), diffs


def substrate_test(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Exploratory. Family centroids (luciferases only) in standardized composition space;
    statistic = mean distance between the coelenterazine families; exact null over every
    choice of the same number of families."""
    luc = table[table.role == "luciferase"]
    cols = ["gravy", "isoelectric_point", "aromaticity"] + [c for c in table if c.startswith("aa_")]
    z = (luc[cols] - luc[cols].mean()) / luc[cols].std()
    fams = [f for f in FAMILIES if f in set(luc.family)]
    centroids = z.groupby(luc["family"]).mean().loc[fams]
    dist = pd.DataFrame(np.linalg.norm(centroids.values[:, None] - centroids.values[None], axis=2),
                        index=fams, columns=fams)
    coel = [f for f in fams if LUCIFERIN[f] == "coelenterazine"]

    def mean_within(group):
        return np.mean([dist.loc[a, b] for a, b in itertools.combinations(group, 2)])

    observed = mean_within(coel)
    null = np.array([mean_within(g) for g in itertools.combinations(fams, len(coel))])
    p = np.mean(null <= observed)
    summary = pd.DataFrame([{"coelenterazine_families": ", ".join(coel), "n_families": len(fams),
                             "observed_mean_distance": observed, "null_median_distance": np.median(null),
                             "label_sets_tested": len(null), "p_exact_one_sided": p}])
    return summary, dist, p


def plot_replication(diffs: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, (feature, direction) in zip(axes, PREDICTED.items()):
        d = diffs[diffs.feature == feature]
        colours = ["#e67e22" if s.startswith("new") else "#1f4e79" for s in d["set"]]
        ax.barh([NAMES[f] for f in d.family], d.difference_sd, color=colours)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(f"{'GRAVY' if feature == 'gravy' else 'Aromaticity'}: luciferase − control "
                     f"(predicted {'lower' if direction < 0 else 'higher'})", fontsize=9.5)
        ax.set_xlabel("difference in SD units", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#1f4e79", "#e67e22")]
    fig.legend(handles, ["original families (hint)", "new families (pre-registered test)"],
               loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(RESULTS_DIR / "figures" / "h2_replication.png", dpi=200)
    plt.close(fig)


def plot_substrate(dist: pd.DataFrame, p: float) -> None:
    order = sorted(dist.index, key=lambda f: (LUCIFERIN[f] != "coelenterazine", FAMILIES.index(f)))
    d = dist.loc[order, order]
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(d.values, cmap="viridis_r")
    labels = [f"{NAMES[f]}{' *' if LUCIFERIN[f] == 'coelenterazine' else ''}" for f in order]
    ax.set_xticks(range(len(order)), labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)), labels, fontsize=8)
    n_coel = sum(LUCIFERIN[f] == "coelenterazine" for f in order)
    ax.add_patch(plt.Rectangle((-0.5, -0.5), n_coel, n_coel, fill=False, ec="red", lw=1.5))
    fig.colorbar(im, ax=ax, shrink=0.8, label="Composition distance (lower = more similar)")
    ax.set_title(f"Family composition distances (* coelenterazine users; exact p = {p:.2f})", fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "substrate_distance.png", dpi=200)
    plt.close(fig)


def run() -> pd.DataFrame:
    table = load()
    tests, diffs = run_tests(table)
    tests.to_csv(RESULTS_DIR / "convergence_tests.csv", index=False, float_format="%.4g", lineterminator="\n")
    diffs.to_csv(RESULTS_DIR / "convergence_by_family.csv", index=False, float_format="%.4f", lineterminator="\n")
    substrate, dist, p = substrate_test(table)
    substrate.to_csv(RESULTS_DIR / "substrate_convergence.csv", index=False, float_format="%.4g", lineterminator="\n")
    plot_replication(diffs)
    plot_substrate(dist, p)
    new = diffs[diffs.family.isin(NEW) & diffs.feature.isin(PREDICTED)]
    passed = (new.observed_direction == new.predicted_direction).all()
    print(new[["family", "feature", "difference", "predicted_direction", "observed_direction"]].to_string(index=False))
    print(f"\nPre-registered replication: {'PASSED' if passed else 'FAILED'}\n")
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(tests.drop(columns=["permutation_statistic"]).to_string(index=False))
    print("\n", substrate.to_string(index=False))
    return tests


if __name__ == "__main__":
    run()

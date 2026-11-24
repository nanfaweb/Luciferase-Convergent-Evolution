"""Novelty phase, part 2: where in the 3D structure does the aromatic excess sit?

For each protein with an AlphaFold model, each residue gets a burial score: the number of
other C-alpha atoms within 10 A. Only confidently predicted residues (pLDDT >= 70) are
used. Within each protein, residues above the protein's median burial count as "buried"
(core and pocket-lining); the rest count as "exposed". This is a simple,
dependency-free proxy for pocket detection.

For buried and exposed residues separately, it computes aromatic % (F, W, Y) and mean
Kyte-Doolittle hydropathy. It then repeats the family-aware luciferase-vs-control tests
from src.convergence. These tests are exploratory (two-sided), not pre-registered.

Model files are cached in data/structures/ (gitignored). Writes results/burial_by_protein.csv,
results/aromatic_burial.csv and results/figures/aromatic_burial.png.

Usage:
    python -m src.pockets   (needs src.features output)
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParamData import kd

from src.align import load_proteins
from src.convergence import NAMES, NEW, permutation_p, with_controls
from src.features import RESULTS_DIR
from src.structure import AFDB_API

STRUCT_DIR = RESULTS_DIR.parent / "data" / "structures"
RADIUS, MIN_PLDDT = 10.0, 70.0
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
                "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
                "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
METRICS = ["aromatic_pct_buried", "aromatic_pct_exposed", "hydropathy_buried", "hydropathy_exposed"]


def model_path(session: requests.Session, accession: str):
    path = STRUCT_DIR / f"{accession}.pdb"
    if path.exists():
        return path
    response = session.get(AFDB_API.format(accession=accession), timeout=60)
    if not (response.ok and response.json()):
        return None
    pdb = session.get(response.json()[0]["pdbUrl"], timeout=120)
    pdb.raise_for_status()
    path.write_bytes(pdb.content)
    return path


def burial_metrics(path) -> dict:
    chain = next(PDBParser(QUIET=True).get_structure("m", path)[0].get_chains())
    residues = [r for r in chain if "CA" in r and r.get_resname() in THREE_TO_ONE]
    coords = np.array([r["CA"].coord for r in residues])
    plddt = np.array([r["CA"].bfactor for r in residues])
    aa = np.array([THREE_TO_ONE[r.get_resname()] for r in residues])
    dist = np.linalg.norm(coords[:, None] - coords[None], axis=2)
    neighbours = (dist < RADIUS).sum(axis=1) - 1
    confident = plddt >= MIN_PLDDT
    cut = np.median(neighbours[confident])
    out = {"residues_confident": int(confident.sum()), "burial_median_cut": cut}
    for name, mask in (("buried", confident & (neighbours > cut)), ("exposed", confident & (neighbours <= cut))):
        out[f"aromatic_pct_{name}"] = 100 * np.isin(aa[mask], list("FWY")).mean()
        out[f"hydropathy_{name}"] = np.mean([kd[a] for a in aa[mask]])
    return out


def run() -> pd.DataFrame:
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    proteins = load_proteins()
    session = requests.Session()
    rows = []
    for prot in proteins[proteins["source"] == "uniprot"].itertuples():
        path = model_path(session, prot.accession)
        if path is not None:
            rows.append({"accession": prot.accession, "family": str(prot.family), "role": prot.role,
                         **burial_metrics(path)})
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS_DIR / "burial_by_protein.csv", index=False, float_format="%.3f", lineterminator="\n")

    families = with_controls(table)
    diff_rows = []
    for fam in families:
        g = table[table.family == fam]
        row = {"family": fam, "set": "new" if fam in NEW else "original",
               "proteins": f"{(g.role == 'luciferase').sum()} luc / {(g.role == 'control').sum()} ctrl"}
        for m in METRICS:
            row[m] = g.loc[g.role == "luciferase", m].median() - g.loc[g.role == "control", m].median()
        diff_rows.append(row)
    for m in METRICS:
        stat, p = permutation_p(table, m, families, None)
        diff_rows.append({"family": "ALL (permutation, two-sided)", "set": "", "proteins": "",
                          m: p, "statistic_" + m: stat})
    diffs = pd.DataFrame(diff_rows)
    diffs.to_csv(RESULTS_DIR / "aromatic_burial.csv", index=False, float_format="%.4g", lineterminator="\n")

    per_fam = diffs[diffs.family.isin(families)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    y = np.arange(len(per_fam))
    for ax, (kind, label) in zip(axes, [("aromatic_pct", "Aromatic % (luciferase − control)"),
                                        ("hydropathy", "Hydropathy (luciferase − control)")]):
        ax.barh(y - 0.2, per_fam[f"{kind}_buried"], 0.4, color="#1f4e79", label="buried (core/pocket)")
        ax.barh(y + 0.2, per_fam[f"{kind}_exposed"], 0.4, color="#85c1e9", label="exposed (surface)")
        ax.set_yticks(y, [NAMES[f] for f in per_fam.family], fontsize=8)
        ax.axvline(0, color="black", lw=0.8)
        ax.invert_yaxis()
        ax.set_title(label, fontsize=9.5)
        ax.grid(axis="x", alpha=0.3)
        ax.tick_params(labelsize=8)
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Where the luciferase-vs-control difference sits in AlphaFold models (pLDDT ≥ 70)", fontsize=10)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "figures" / "aromatic_burial.png", dpi=200)
    plt.close(fig)
    print(f"AlphaFold models analysed: {len(table)}")
    print(diffs.drop(columns=[c for c in diffs if c.startswith("statistic_")]).to_string(index=False))
    return diffs


if __name__ == "__main__":
    run()

"""Build the final report as a PDF, with figures and tables taken from results/.

Run after the analysis modules. Writes three summary figures to
results/figures/ and the report to docs/FINAL_REPORT.pdf.

Usage:
    python -m src.report_pdf
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, Preformatted,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

from src.download import FAMILIES
from src.features import RESULTS_DIR

ROOT = RESULTS_DIR.parent
FIG = RESULTS_DIR / "figures"
OUT = ROOT / "docs" / "FINAL_REPORT.pdf"
NAMES = {"firefly": "Firefly", "renilla_type": "Renilla-type", "bacterial_luxa": "Bacterial LuxA",
         "fungal_luz": "Fungal Luz", "copepod": "Copepod", "dinoflagellate": "Dinoflagellate",
         "ostracod": "Ostracod", "click_beetle": "Click beetle"}
ACCENT = colors.HexColor("#1f4e79")

# ---------------------------------------------------------------- fonts / styles
FONT_DIR = Path("C:/Windows/Fonts")
if (FONT_DIR / "arial.ttf").exists():
    for name, file in [("Body", "arial.ttf"), ("Body-Bold", "arialbd.ttf"),
                       ("Body-Italic", "ariali.ttf"), ("Body-BoldItalic", "arialbi.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(FONT_DIR / file)))
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold",
                                  italic="Body-Italic", boldItalic="Body-BoldItalic")
    BODY, BOLD = "Body", "Body-Bold"
else:  # fallback: built-in fonts (no Unicode minus etc.)
    BODY, BOLD = "Helvetica", "Helvetica-Bold"

ss = getSampleStyleSheet()
S = {
    "title": ParagraphStyle("t", parent=ss["Title"], fontName=BOLD, fontSize=19, leading=24, textColor=ACCENT),
    "sub": ParagraphStyle("s", parent=ss["Normal"], fontName=BODY, fontSize=11, alignment=TA_CENTER, leading=15),
    "h1": ParagraphStyle("h1", parent=ss["Heading1"], fontName=BOLD, fontSize=15, textColor=ACCENT, spaceBefore=10, spaceAfter=6),
    "h2": ParagraphStyle("h2", parent=ss["Heading2"], fontName=BOLD, fontSize=12, textColor=ACCENT, spaceBefore=8, spaceAfter=4),
    "h3": ParagraphStyle("h3", parent=ss["Heading3"], fontName=BOLD, fontSize=10.5, spaceBefore=6, spaceAfter=2),
    "p": ParagraphStyle("p", parent=ss["Normal"], fontName=BODY, fontSize=9.8, leading=13.5, spaceAfter=5),
    "cell": ParagraphStyle("c", parent=ss["Normal"], fontName=BODY, fontSize=8.3, leading=10.3),
    "cellb": ParagraphStyle("cb", parent=ss["Normal"], fontName=BOLD, fontSize=8.3, leading=10.3, textColor=colors.white),
    "cap": ParagraphStyle("cap", parent=ss["Normal"], fontName=BODY, fontSize=8.5, leading=11,
                          textColor=colors.HexColor("#444444"), spaceAfter=8, alignment=TA_CENTER),
    "ref": ParagraphStyle("r", parent=ss["Normal"], fontName=BODY, fontSize=8.6, leading=11.2, spaceAfter=3,
                          leftIndent=14, firstLineIndent=-14),
    "box": ParagraphStyle("box", parent=ss["Normal"], fontName=BODY, fontSize=8.6, leading=10.8, alignment=TA_CENTER),
}


def P(text: str, style: str = "p") -> Paragraph:
    return Paragraph(text, S[style])


def bullets(items: list[str]) -> list:
    return [Paragraph(i, S["p"], bulletText="•") for i in items]


def table(rows: list[list], widths: list[float], bold_last: bool = False) -> Table:
    data = [[Paragraph(str(c), S["cellb"]) for c in rows[0]]]
    data += [[Paragraph(str(c), S["cell"]) for c in r] for r in rows[1:]]
    t = Table(data, colWidths=[w * cm for w in widths], repeatRows=1)
    style = [("BACKGROUND", (0, 0), (-1, 0), ACCENT),
             ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef3f8")]),
             ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#b0b8c4")),
             ("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
    if bold_last:
        style.append(("BACKGROUND", (0, -1), (-1, -1), colors.HexColor("#d6e2ef")))
    t.setStyle(TableStyle(style))
    return t


def figure(path: Path, width_cm: float, caption: str) -> KeepTogether:
    img = Image(str(path))
    img.drawWidth, img.drawHeight = width_cm * cm, width_cm * cm * img.imageHeight / img.imageWidth
    return KeepTogether([img, P(caption, "cap")])


def tabcap(text: str) -> Paragraph:
    return P(f"<b>{text}</b>", "cap")


def captioned(caption: str, rows: list[list], widths: list[float], **kw) -> list:
    return [KeepTogether([tabcap(caption), table(rows, widths, **kw)])]


# ---------------------------------------------------------------- data
def load() -> dict:
    d = {name: pd.read_csv(RESULTS_DIR / f"{name}.csv", keep_default_na=False) for name in
         ["pairwise", "control_comparison", "stats", "blast_summary", "motif_hits", "prosite_patterns",
          "go_by_family", "tryptic_summary", "features_by_family"]}
    d["alphafold"] = pd.read_csv(RESULTS_DIR / "alphafold_confidence.csv")
    d["peptides"] = pd.read_csv(RESULTS_DIR / "tryptic_peptides.csv")
    d["blast"] = pd.read_csv(RESULTS_DIR / "blast_all_vs_all.tsv", sep="\t")
    d["manifest"] = pd.read_csv(ROOT / "data" / "curated" / "manifest.csv", keep_default_na=False)
    pw = d["pairwise"].copy()
    beetle = pw["family_a"].isin(["firefly", "click_beetle"]) & pw["family_b"].isin(["firefly", "click_beetle"])
    pw.loc[(pw["category"] == "cross-family") & beetle, "category"] = "firefly vs click beetle"
    d["pairs"] = pw
    return d


ORDER = ["same family", "luciferase vs own control", "firefly vs click beetle", "cross-family"]
LABEL = {"same family": "Same family", "luciferase vs own control": "Luciferase vs own control",
         "firefly vs click beetle": "Firefly vs click beetle", "cross-family": "Other cross-family"}


def make_figures(d: dict) -> None:
    pw = d["pairs"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, col, title in [(axes[0], "global_identity_pct", "Global identity (%)"),
                           (axes[1], "local_zscore", "Smith–Waterman z-score")]:
        data = [pw.loc[pw["category"] == c, col] for c in ORDER]
        ax.boxplot(data, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="#9ecae1", color="#1f4e79"), medianprops=dict(color="#c0392b"))
        ax.set_xticks(range(1, 5), [LABEL[c].replace(" vs ", "\nvs ") for c in ORDER], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    axes[1].axhline(3, ls="--", color="grey", lw=0.8)
    axes[1].text(4.45, 3, "z = 3", fontsize=7, va="bottom", ha="right", color="grey")
    fig.tight_layout()
    fig.savefig(FIG / "similarity_by_category.png", dpi=200)
    plt.close(fig)

    cc = d["control_comparison"]
    feats = [("gravy", "GRAVY (luciferase − control)"), ("aromaticity", "Aromaticity (luciferase − control)")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.2))
    for ax, (feat, title) in zip(axes, feats):
        sub = cc[cc["feature"] == feat]
        vals = sub["difference"].astype(float)
        ax.barh([NAMES[f] for f in sub["family"]], vals,
                color=["#c0392b" if v < 0 else "#2e86c1" for v in vals])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(title, fontsize=10)
        ax.tick_params(labelsize=8)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "control_differences.png", dpi=200)
    plt.close(fig)

    af = d["alphafold"].dropna(subset=["mean_plddt"])
    fams = [f for f in FAMILIES if f in set(af["family"])]
    fig, ax = plt.subplots(figsize=(10, 3.4))
    for i, fam in enumerate(fams):
        for role, colour, dx in [("luciferase", "#1f4e79", -0.12), ("control", "#e67e22", 0.12)]:
            v = af[(af["family"] == fam) & (af["role"] == role)]["mean_plddt"]
            ax.scatter([i + dx] * len(v), v, s=18, color=colour, alpha=0.8,
                       label=role if i == 0 else None)
    ax.axhspan(90, 100, color="#2e86c1", alpha=0.08)
    ax.axhspan(70, 90, color="#85c1e9", alpha=0.08)
    ax.set_xticks(range(len(fams)), [NAMES[f] for f in fams], fontsize=8)
    ax.set_ylabel("Mean pLDDT", fontsize=9)
    ax.set_ylim(30, 100)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("AlphaFold confidence per protein (>90 very high, 70–90 confident)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "alphafold_plddt.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- report
def build(d: dict) -> list:
    m, pw = d["manifest"], d["pairs"]
    n_luc, n_ctl = (m["role"] == "luciferase").sum(), (m["role"] == "control").sum()
    cat = pw.groupby("category").agg(n=("global_identity_pct", "size"),
                                     ident=("global_identity_pct", "median"), z=("local_zscore", "median"))
    other = pw[pw["category"] == "cross-family"]
    n_z3, max_z = int((other["local_zscore"] > 3).sum()), other["local_zscore"].max()
    blast = d["blast"]
    cross = blast[blast["query_family"] != blast["subject_family"]]
    cross_fams = sorted({tuple(sorted(x)) for x in zip(cross["query_family"], cross["subject_family"])})
    pep = d["peptides"]
    shared_peps = pep[pep["n_families_with_peptide"] > 1]["peptide"].nunique()
    ari = float(d["stats"].iloc[0]["statistic"])
    sign = d["stats"][d["stats"]["test"].str.startswith("Sign test")]
    kw = d["stats"][d["stats"]["test"].str.startswith("Kruskal")]
    af = d["alphafold"]
    s: list = []

    # title page
    s += [Spacer(1, 3 * cm),
          P("Molecular Signatures of Convergent Evolution", "title"),
          P("Comparative Sequence Analysis of Independently Evolved Luciferases", "sub"),
          Spacer(1, 0.6 * cm),
          P("Final Report · Bioinformatics course project · FAST-NU Lahore", "sub"),
          P(date.today().strftime("%B %Y"), "sub"),
          Spacer(1, 1.2 * cm),
          P("<b>Abstract.</b> Bioluminescence evolved independently many times, and each time a different "
            "enzyme, a luciferase, was recruited to produce light. We asked whether these independently "
            f"evolved luciferases share sequence-level or physicochemical features. We curated {len(m)} "
            f"proteins: {n_luc} luciferases from {m.loc[m.role == 'luciferase', 'family'].nunique()} families "
            f"(9 independent origins of light) and {n_ctl} non-luminous relatives used as controls. We then "
            "applied the course's methods: database curation, Gene Ontology, global and local alignment by "
            "dynamic programming, dot plots, BLAST, phylogenetic trees, motif search, an in-silico trypsin "
            "digest, physicochemical profiling, clustering, statistics and AlphaFold structure confidence. "
            "Independent luciferase families share <b>no</b> detectable sequence similarity, motif or peptide; "
            "the only cross-family similarity (firefly with click beetle) reflects a shared ancestral enzyme "
            "family. Physicochemical profiles are dominated by family ancestry (adjusted Rand index "
            f"{ari:.2f}). In all 5 families with controls, luciferases are less hydrophobic and more aromatic "
            "than their non-luminous relatives (sign test p = 0.0625): a consistent hint that is not "
            "statistically significant. Luciferases converged in function, not in sequence."),
          PageBreak()]

    # 1 problem
    s += [P("1. Problem formulation", "h1"),
          P("<b>Background.</b> Bioluminescence, the production of light by living organisms, has evolved "
            "independently many times: in beetles, marine animals, bacteria, fungi, crustaceans and "
            "dinoflagellates. In each case an enzyme called a <b>luciferase</b> drives a light-emitting "
            "oxidation reaction. Different lineages recruited luciferases from unrelated protein families: "
            "acyl-CoA synthetases (ACS) in beetles, haloalkane dehalogenases (HLD) in the sea pansy "
            "<i>Renilla</i>, and flavin-dependent monooxygenases in bacteria. These enzymes perform the same "
            "<i>function</i>. Whether they also converged at the <i>molecular</i> level is an open, testable "
            "question."),
          P("<b>Research question.</b> Do independently evolved luciferases share common protein-sequence or "
            "physicochemical characteristics despite their different evolutionary origins?"),
          P("<b>Hypotheses.</b>")]
    s += bullets(["<b>H1 (sequence).</b> Luciferases from independent origins share no detectable sequence "
                  "similarity, motifs or peptides.",
                  "<b>H2 (properties).</b> If light production imposes shared constraints, luciferases should "
                  "differ from their closest non-luminous relatives (controls) <b>in the same direction</b> "
                  "across independent families.",
                  "<b>H0 for H2.</b> Properties are set by protein-family ancestry: luciferases group by family, "
                  "and luciferase-vs-control differences point in random directions."])
    s += [P("<b>Scope.</b> A protein-level computational comparison of public sequences. No DNA motif search, "
            "no wet-lab work, and no claim that any property <i>causes</i> light emission.")]

    # 2 existing
    s += [P("2. Existing solutions and related work", "h1"), P("2.1 Studies of luciferase evolution", "h2")]
    s += bullets(["Delroisse et al. [1] compared luciferases and photoproteins across many phyla and described "
                  "several independent evolutionary groups.",
                  "He et al. [2] showed that beetle bioluminescence has multiple origins: fireflies and click "
                  "beetles each recruited luciferase from the ACS family.",
                  "Lau et al. [4] showed that a brittle star recruited a haloalkane dehalogenase as a luciferase "
                  "independently of <i>Renilla</i>, a case of parallel evolution.",
                  "Kotlobay et al. [3] identified the fungal luciferase Luz as a distinct protein family; "
                  "Takenaka et al. [5] and Liu et al. [6] characterised copepod and dinoflagellate luciferases."])
    s += [P("2.2 Standard computational methods", "h2"),
          table([["Task", "Standard solution", "Used here"],
                 ["Pairwise alignment", "Needleman–Wunsch [7], Smith–Waterman [8], BLOSUM62 [9]", "Yes (E2, E3)"],
                 ["Homology search", "BLAST+ [10]", "Yes (E5)"],
                 ["Distance trees", "Neighbour joining [11]", "Yes (E6)"],
                 ["Physicochemical profiling", "ExPASy ProtParam [13], Kyte–Doolittle hydropathy [14]", "Yes (E1)"],
                 ["Motif databases", "PROSITE [18]", "Yes (E7)"],
                 ["Domain annotation", "Profile HMMs via Pfam/InterPro", "For curation"],
                 ["Structure prediction", "AlphaFold [15], AlphaFold DB [16]", "Yes (E10)"]],
                [4, 8.5, 3.5]),
          Spacer(1, 6),
          P("<b>Gap addressed.</b> Earlier studies mostly built family trees or characterised single enzymes. "
            "We combine the standard methods on <b>one curated, multi-origin dataset with matched non-luminous "
            "controls</b>, which separates traits inherited from the protein family from traits that could be "
            "shared because the protein makes light.")]

    # 3 methodology
    s += [PageBreak(), P("3. Proposed solution: methodology and framework", "h1"),
          P("The framework is a pipeline from public databases to family-level statistics (Figure 1). Every "
            "step is a Python module in <font face='Courier'>src/</font> and the whole pipeline re-runs from "
            "the notebook <font face='Courier'>notebooks/01_analysis.ipynb</font>.")]
    box = lambda t: Paragraph(t, S["box"])
    arrow = Paragraph("▼", S["box"])
    stages = [
        "<b>1. Data collection</b><br/>UniProt and NCBI Protein REST APIs",
        f"<b>2. Curation</b> (rules in seeds.csv)<br/>{len(m)} proteins: {n_luc} luciferases in 8 families + {n_ctl} controls",
    ]
    branches = [["<b>Sequence</b><br/>NW / SW alignment,<br/>shuffle z-scores,<br/>dot plots, NJ trees",
                 "<b>Homology</b><br/>BLAST+<br/>all-vs-all",
                 "<b>Motifs and peptides</b><br/>PROSITE multi-pattern<br/>scan, trypsin digest",
                 "<b>Properties</b><br/>ProtParam: GRAVY,<br/>pI, aromaticity,<br/>composition",
                 "<b>Structure</b><br/>AlphaFold DB<br/>confidence"]]
    flow = [[box(stages[0])], [arrow], [box(stages[1])], [arrow]]
    t1 = Table(flow, colWidths=[16 * cm])
    t1.setStyle(TableStyle([("BOX", (0, 0), (0, 0), 0.8, ACCENT), ("BOX", (0, 2), (0, 2), 0.8, ACCENT),
                            ("BACKGROUND", (0, 0), (0, 0), colors.HexColor("#eef3f8")),
                            ("BACKGROUND", (0, 2), (0, 2), colors.HexColor("#eef3f8"))]))
    t2 = Table([[box(b) for b in branches[0]]], colWidths=[3.2 * cm] * 5)
    t2.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.8, ACCENT), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f7f9fc"))]))
    t3 = Table([[arrow], [box("<b>3. Comparison and statistics</b><br/>PCA · Ward clustering vs family (ARI) · "
                              "Kruskal–Wallis · luciferase-vs-own-control sign test (unit = family)")],
                [arrow], [box("<b>4. Answer</b>: shared sequence? shared property shift?")]],
               colWidths=[16 * cm])
    t3.setStyle(TableStyle([("BOX", (0, 1), (0, 1), 0.8, ACCENT), ("BOX", (0, 3), (0, 3), 0.8, ACCENT),
                            ("BACKGROUND", (0, 1), (0, 1), colors.HexColor("#eef3f8")),
                            ("BACKGROUND", (0, 3), (0, 3), colors.HexColor("#d6e2ef"))]))
    s += [KeepTogether([t1, t2, t3, P("Figure 1. Analysis framework.", "cap")]),
          P("3.1 Key design choices", "h2")]
    s += bullets(["<b>Matched controls.</b> Each of 5 families has 1–3 non-luminous homologs from the same protein "
                  "family. Comparing a luciferase with its own control cancels out family-level traits.",
                  "<b>The family is the unit of evidence.</b> Sequences within a family share ancestry, so the "
                  "statistics compare families, not individual proteins.",
                  "<b>No mixed-family tree.</b> Unrelated proteins cannot be meaningfully aligned together, so trees "
                  "are built within families only.",
                  "<b>Composition, not size.</b> Length and weight reflect domain architecture, so cross-family "
                  "comparisons use composition-type features.",
                  "<b>Chance baselines.</b> Alignment scores are judged against shuffled sequences (z-scores), and "
                  "BLAST hits against an E-value threshold, so every 'similar' claim has a null model."])

    # 4 implementation
    s += [PageBreak(), P("4. Implementation and experiments", "h1"), P("4.1 Environment", "h2")]
    s += bullets(["Python 3.13 with Biopython 1.86 [12], pandas, NumPy, SciPy, Matplotlib and scikit-learn.",
                  "Local BLAST+ 2.17.0 (run in a temporary folder because BLAST cannot handle spaces in paths).",
                  "Data from the UniProt [17], NCBI Protein, PROSITE and AlphaFold DB REST APIs."])
    fam_rows = [["Family", "Luciferases", "Controls", "Origin of the luciferase"]]
    origin = {"firefly": "ACS family", "click_beetle": "ACS family, independent of fireflies [2]",
              "renilla_type": "HLD family, two independent recruitments (Renilla, brittle star) [4]",
              "bacterial_luxa": "Luciferase-like monooxygenase family", "fungal_luz": "Luz family [3]",
              "copepod": "Unique family [5]", "dinoflagellate": "Unique three-domain enzyme [6]",
              "ostracod": "Cypridinid luciferase"}
    for fam in FAMILIES:
        sub = m[m["family"] == fam]
        ctl = sub[sub["role"] == "control"]
        fam_rows.append([NAMES[fam], (sub["role"] == "luciferase").sum(),
                         f"{len(ctl)}" + (f" ({', '.join(ctl['label'])})" if len(ctl) else "") if len(ctl) else "–",
                         origin[fam]])
    fam_rows.append(["<b>Total</b>", f"<b>{n_luc}</b>", f"<b>{n_ctl}</b>", "<b>9 independent origins</b>"])
    s += [P("4.2 Dataset", "h2"),
          P("Selection rules: functional evidence that the protein is a luciferase; complete sequence (6 "
            "near-complete clones kept and flagged); no two sequences in a family ≥97% identical; expert-reviewed "
            "entries preferred; every record versioned, hashed and re-verified against live UniProt."),
          tabcap("Table 1. Curated dataset."),
          table(fam_rows, [2.8, 2, 5.6, 5.6], bold_last=True)]

    s += [P("4.3 Experiments by course topic", "h2"),
          P("Each subsection names the course topic, what was done with it in this project, and the settings. "
            "Experiment numbers (E1–E10) link to the results in Section 5.")]
    topics = [
        ("Biological databases and Gene Ontology (Week 3)",
         "Sequences, metadata and annotations were fetched programmatically from UniProtKB and NCBI Protein. "
         "Each record's evidence level, reviewed status, InterPro/Pfam domains and GO terms were stored in "
         "<font face='Courier'>data/curated/manifest.csv</font>. GO molecular-function terms were then summarised "
         "per family and role (<font face='Courier'>results/go_by_family.csv</font>) to confirm that each "
         "luciferase carries a luciferase/monooxygenase activity while its controls carry the ancestral activity."),
        ("Proteomics: physicochemical properties (Week 14) — E1",
         "Biopython's ProtParam (the ExPASy method [13]) computed length, molecular weight, isoelectric point, "
         "GRAVY hydropathy (Kyte–Doolittle [14]), aromaticity and the percentage of each amino acid for every "
         "protein."),
        ("Similarity and homology (Week 6); global and local alignment (Week 7) — E2",
         "All 2,211 protein pairs were aligned with Biopython's PairwiseAligner using BLOSUM62 [9] and gap "
         "open −10 / extend −0.5. Global (Needleman–Wunsch [7]) alignment gave percent identity over the shorter "
         "sequence. Local (Smith–Waterman [8]) alignment gave a score that was compared with 20 shuffled copies of "
         "the second sequence, producing a z-score: z ≈ 0 means chance; z > 3 is borderline; large z means homology."),
        ("Dynamic programming (Week 8) — E3",
         "A Needleman–Wunsch function was written by hand: it fills the score table cell by cell "
         "(match/mismatch from BLOSUM62, linear gap −8), then traces back from the bottom-right cell to recover "
         "the alignment. It was run on the textbook pair HEAGAWGHEE / PAWHEAE and checked against Biopython."),
        ("Dot plots (Week 9) — E4",
         "Three dot plots compare firefly luciferase with a second firefly, with the non-luminous click beetle "
         "AbLL, and with bacterial LuxA. A dot is drawn where a 10-residue window scores ≥23 with BLOSUM62."),
        ("BLAST and multiple pattern matching (Week 12) — E5, E7",
         "A local BLAST+ database of all 67 proteins was built and searched all-vs-all with blastp (E < 10<super>−3</super>). "
         "For motif search, every PROSITE pattern [18] that UniProt lists for any protein in the dataset was "
         "converted to a regular expression and all patterns were scanned against all sequences at once "
         "(multiple pattern matching)."),
        ("Motif search (Week 2) — E7",
         "The PROSITE scan above is the protein-level form of motif finding: it asks whether any short, conserved "
         "signature appears in more than one independent luciferase family."),
        ("Phylogenetics (Week 10) — E6",
         "Within each family, global identities were turned into Poisson distances (−ln(1−p)) and a "
         "neighbour-joining tree [11] was built, rooted on the control least similar to the luciferases (or at the "
         "midpoint where there is no control)."),
        ("Antibiotics sequencing / mass spectrometry (Week 5) and peptide sequencing (Week 14) — E8",
         "Each protein was digested in silico with trypsin (cut after K or R unless the next residue is P). "
         "Peptides of 6–30 residues, the range a mass spectrometer usually detects, were given a monoisotopic "
         "mass. Shared peptides would reveal shared sequence; unique peptides form a mass fingerprint per family."),
        ("Clustering (Week 11) and statistics — E9",
         "Composition features were standardised and reduced by PCA. Ward hierarchical clustering into 8 clusters "
         "was compared with the true families by the adjusted Rand index (ARI [20]). Kruskal–Wallis tests [19] "
         "compared the 8 families. For H2, the luciferase median was compared with its own control median in each "
         "of the 5 families, and a two-sided sign test asked whether the direction was consistent."),
        ("Hidden Markov Models (Week 13)",
         "Pfam and InterPro domain assignments, which are built from profile HMMs, were used during curation to "
         "confirm the family of every protein. No HMM was trained, because unrelated luciferases cannot be "
         "aligned into one profile."),
        ("Deep learning for biomolecules (Week 15) — E10",
         "For every UniProt accession, the AlphaFold DB API [16] returned the predicted model's mean pLDDT "
         "(AlphaFold's per-residue confidence [15], 0–100), giving the structural coverage of the dataset."),
    ]
    for i, (head, body) in enumerate(topics, 1):
        s += [P(f"4.3.{i} {head}", "h3"), P(body)]

    # 5 results
    s += [PageBreak(), P("5. Results", "h1"), P("5.1 Sequence similarity (E2, E4, E5)", "h2"),
          P("This is the headline result. Pairs within a family are clearly homologous; pairs between independent "
            "families score at chance level (Table 2, Figure 2).")]
    rows = [["Pair type", "Pairs", "Median identity", "Median SW z-score"]]
    for c in ORDER:
        rows.append([LABEL[c], f"{int(cat.loc[c, 'n']):,}", f"{cat.loc[c, 'ident']:.0f}%", f"{cat.loc[c, 'z']:.1f}"])
    s += [*captioned("Table 2. Pairwise alignment by pair type.", rows, [6, 2.5, 3.5, 4]), Spacer(1, 8),
          figure(FIG / "similarity_by_category.png", 16,
                 "Figure 2. Identity and Smith–Waterman z-score by pair type (dashed line: z = 3).")]
    s += bullets([f"<b>BLAST.</b> {len(blast)} significant hits. The only cross-family hits link "
                  + ", ".join(" and ".join(NAMES[x] for x in pair) for pair in cross_fams)
                  + " (best E-value below 10<super>−126</super>, reported as 0), which are expected ACS relatives. No other "
                  "luciferase found a luciferase from another family.",
                  f"<b>Weak alignment hits.</b> {n_z3} of {len(other):,} unrelated pairs reached z > 3 (maximum "
                  f"{max_z:.1f}). With 20 shuffles and over 1,000 tests a few are expected by chance, and BLAST "
                  "supported none of them.",
                  "<b>Dot plots.</b> Same family: a continuous diagonal. Luciferase vs control: a broken diagonal. "
                  "Firefly vs LuxA: scattered noise only (Figure 3)."])
    s += [figure(FIG / "identity_heatmap.png", 13.5,
                 "Figure 3a. Pairwise global identity for all 67 proteins, ordered by family (white lines). "
                 "Bright blocks appear only on the diagonal, plus the firefly/click-beetle block."),
          figure(FIG / "dotplots.png", 16, "Figure 3b. Dot plots: same family, luciferase vs control, unrelated families.")]

    dp = (RESULTS_DIR / "dp_example.txt").read_text(encoding="utf-8").strip()
    s += [KeepTogether([P("5.2 Dynamic programming worked example (E3)", "h2"),
          P("The hand-written table reaches the same optimal score as Biopython for identical settings, "
            "confirming the implementation."),
          Preformatted(dp, ParagraphStyle("pre", fontName="Courier", fontSize=7.6, leading=9.2,
                                          backColor=colors.HexColor("#f4f6f8"), borderPadding=5,
                                          leftIndent=5, spaceAfter=8))])]

    s += [P("5.3 Phylogenetic trees (E6)", "h2"),
          P("The trees reproduce known biology: firefly subfamilies, <i>Vibrio</i> vs <i>Photobacterium</i> LuxA, "
            "and the <i>Armillaria</i> Luz cluster. In the <i>Renilla</i>-type tree the brittle-star luciferase "
            "groups with the <b>brittle star's own non-luminous dehalogenase</b>, not with <i>Renilla</i> "
            "luciferase, directly showing two separate recruitments of the same enzyme family."),
          figure(FIG / "trees.png", 16, "Figure 4. Within-family neighbour-joining trees.")]

    s += [P("5.4 Motifs and Gene Ontology (E7, Week 3)", "h2")]
    mh = d["motif_hits"]
    rows = [["PROSITE pattern", "Description", "Families hit (luciferase / control)", "Hits"]]
    for pid, g in mh.groupby("prosite_id"):
        pat = d["prosite_patterns"].set_index("prosite_id").loc[pid]
        hits = "; ".join(f"{NAMES[f]}: {(gg.role == 'luciferase').sum()} / {(gg.role == 'control').sum()}"
                         for f, gg in g.groupby("family"))
        rows.append([f"{pid}<br/>{pat['name']}", pat["description"], hits, len(g)])
    s += [*captioned("Table 3. PROSITE motif hits.", rows, [4.8, 4.2, 5.4, 1.6]), Spacer(1, 6)]
    s += bullets(["The <b>AMP-binding motif</b> occurs in firefly and click-beetle luciferases <b>and</b> in their "
                  "non-luminous controls: it is inherited from the ACS family, not linked to light.",
                  "The <b>bacterial-luciferase signature</b> occurs in LuxA but not in the SsuD controls.",
                  "<b>No motif is shared between independent luciferase families.</b>"])
    go = d["go_by_family"]
    rows = [["Family", "Role", "GO molecular function (summary)"]]
    for r in go.itertuples():
        rows.append([NAMES.get(r.family, r.family), r.role, r.go_molecular_function or "–"])
    s += [*captioned("Table 4. GO molecular-function terms by family.", rows, [3, 2, 11])]

    s += [P("5.5 Tryptic peptides (E8)", "h2")]
    ts = d["tryptic_summary"]
    rows = [["Family", "Role", "Proteins", "Peptides (6–30 aa)", "Shared with another family"]]
    for r in ts.itertuples():
        rows.append([NAMES[r.family], r.role, r.proteins, r.peptides, r.shared_with_other_family])
    s += [P(f"Of {len(pep):,} tryptic peptides, only {shared_peps} distinct peptides occur in more than one family, "
            "and they are shared only by fireflies and click beetles (Table 5). Every independent family has its own "
            "peptide fingerprint, which in principle could identify it by mass spectrometry."),
          *captioned("Table 5. In-silico trypsin digest.", rows, [3.5, 2.5, 2.5, 3.5, 4])]

    s += [P("5.6 Physicochemical properties, clustering and statistics (E1, E9)", "h2")]
    fb = d["features_by_family"]
    rows = [["Family", "Role", "n", "Length", "pI", "GRAVY", "Aromaticity"]]
    for r in fb.itertuples():
        rows.append([NAMES[r.family], r.role, r.n, f"{float(r.length_median):.0f}",
                     f"{float(r.isoelectric_point_median):.2f}", f"{float(r.gravy_median):.3f}",
                     f"{float(r.aromaticity_median):.3f}"])
    s += [*captioned("Table 6. Median properties by family and role.", rows, [3, 2.2, 1, 2, 2, 2.4, 2.4]),
          Spacer(1, 6)]
    s += bullets([f"<b>Clustering.</b> Ward clustering of the composition features recovers the families with "
                  f"ARI = {ari:.2f} (1 is perfect). Most of the remaining mismatch is fireflies and click beetles "
                  "mixing, as expected for relatives (Figure 5).",
                  "<b>Family differences.</b> Every main property differs strongly between families (Kruskal–Wallis "
                  f"p ≤ {kw['p_value'].astype(float).max():.0e})."])
    s += [figure(FIG / "pca_features.png", 12, "Figure 5. PCA of composition features, coloured by family."),
          figure(FIG / "feature_boxplots.png", 16, "Figure 6. Property distributions by family.")]
    feat_name = {"gravy": "GRAVY (hydrophobicity)", "aromaticity": "Aromaticity", "aa_C_pct": "Cysteine %",
                 "isoelectric_point": "pI"}
    rows = [["Property", "Luciferase vs own control, across 5 families", "Sign test p"]]
    for r in sign.itertuples():
        hi = int(float(r.statistic))
        rows.append([feat_name.get(r.feature, r.feature), f"higher in {hi}, lower in {int(r.n) - hi}", r.p_value])
    s += [P("<b>Luciferase vs own control (H2).</b> In all 5 families with controls, luciferases are less "
            "hydrophobic and more aromatic than their non-luminous relatives (Table 7, Figure 7). With 5 families, "
            "p = 0.0625 is the smallest value the sign test can give, so the result is consistent but not "
            "significant."),
          *captioned("Table 7. Direction of luciferase-vs-control differences.", rows, [5, 7, 3]), Spacer(1, 6),
          figure(FIG / "control_differences.png", 16,
                 "Figure 7. Luciferase minus control median, per family (red: lower, blue: higher).")]

    found = af.dropna(subset=["mean_plddt"])
    rows = [["Family", "Models", "Median pLDDT (luciferases)", "Lowest"]]
    for fam in FAMILIES:
        v = found[(found["family"] == fam) & (found["role"] == "luciferase")]["mean_plddt"]
        if len(v):
            rows.append([NAMES[fam], len(v), f"{v.median():.1f}", f"{v.min():.1f}"])
    s += [P("5.7 AlphaFold structure confidence (E10)", "h2"),
          P(f"{len(found)} of {len(af)} proteins have AlphaFold models; the other {len(af) - len(found)} are the "
            "NCBI-only sequences. Most families are predicted with very high confidence. Copepod luciferases, small "
            "secreted and disulfide-rich proteins, are predicted much less confidently (Table 8, Figure 8)."),
          *captioned("Table 8. AlphaFold confidence by family.", rows, [4, 2.5, 5, 3]), Spacer(1, 6),
          figure(FIG / "alphafold_plddt.png", 16, "Figure 8. Mean pLDDT per protein, by family and role.")]

    # 6 conclusion
    s += [PageBreak(), P("6. Conclusion", "h1")]
    s += bullets(["<b>H1 is supported: no shared sequence.</b> Independent luciferase families share no detectable "
                  "sequence similarity, motifs or peptides. Alignment, dot plots, BLAST, PROSITE and peptide analysis "
                  "agree. The only cross-family similarity (firefly with click beetle) reflects a shared ancestral "
                  "enzyme family, not convergence.",
                  f"<b>Properties are dominated by ancestry.</b> Composition sorts proteins by family (ARI {ari:.2f}), "
                  "and the families differ strongly from each other.",
                  "<b>H2 has a consistent but unproven hint.</b> In all 5 testable families luciferases are less "
                  "hydrophobic and more aromatic than their non-luminous relatives (p = 0.0625 each). It does not "
                  "survive correction for testing 4 properties, and fireflies and click beetles are not fully "
                  "independent, so it is a hypothesis for future work, not a finding.",
                  "<b>Convergence also happens within protein families.</b> Brittle star and <i>Renilla</i>, and "
                  "fireflies and click beetles, each recruited the same enzyme family independently."])
    s += [P("<b>Answer.</b> Independently evolved luciferases converged in <b>function</b>, not in <b>sequence</b>. "
            "Their physicochemical profiles are set mostly by ancestry, with a weak, consistent tendency towards "
            "lower hydrophobicity and higher aromaticity that deserves testing on more families."),
          P("6.1 Limitations", "h2")]
    s += bullets(["Only 5 families have controls, and some groups are small (3 sequences).",
                  "Several sequences are unreviewed or predicted; evidence levels are recorded.",
                  "Secreted luciferases are analysed with their signal peptides.",
                  "Whole-protein features cannot detect active-site or 3D-shape convergence."])
    s += [P("6.2 Future work", "h2")]
    s += bullets(["Compare AlphaFold structures, for example the shape of the active-site pocket.",
                  "Test the hydrophobicity and aromaticity hint on more luciferase families with controls.",
                  "Analyse mature proteins with signal peptides removed."])

    rows = [["Week", "Topic", "Used in this project"],
            ["1", "Intro, origin of replication", "Background only (genome-level; not applicable)"],
            ["2", "Motif search", "PROSITE motif scan (4.3.7, E7)"],
            ["3", "Biological databases, GO", "UniProt/NCBI curation, InterPro/Pfam, GO summary (4.3.1)"],
            ["4", "Genome assembly", "Not applicable (no genome data)"],
            ["5", "Antibiotics sequencing, mass spectrometry", "Trypsin digest and peptide masses (4.3.9, E8)"],
            ["6", "Similarity and homology", "Identity, z-scores and BLAST (4.3.3, E2, E5)"],
            ["7", "Global and local alignment", "Needleman–Wunsch and Smith–Waterman (4.3.3, E2)"],
            ["8", "Dynamic programming", "Hand-written DP table and traceback (4.3.4, E3)"],
            ["9", "Dot plots (sorting by reversals: n/a)", "Dot plots (4.3.5, E4)"],
            ["10", "Phylogenetics", "Neighbour-joining trees (4.3.8, E6)"],
            ["11", "Clustering (transcriptomics: n/a)", "PCA, Ward clustering, ARI (4.3.10, E9)"],
            ["12", "Multiple pattern matching, BLAST", "Multi-pattern PROSITE scan, BLAST+ (4.3.6, E5, E7)"],
            ["13", "Hidden Markov Models", "Pfam/InterPro profile-HMM domains for curation (4.3.11)"],
            ["14", "Proteomics, peptide sequencing", "ProtParam properties, tryptic peptides (4.3.2, 4.3.9)"],
            ["15", "Deep learning for biomolecules", "AlphaFold structure confidence (4.3.12, E10)"]]
    s += [P("6.3 Course coverage", "h2"), tabcap("Table 9. Course topics and where they are used."),
          table(rows, [1.3, 5.7, 9])]

    # references (from the markdown report, the single source)
    md = (ROOT / "docs" / "FINAL_REPORT.md").read_text(encoding="utf-8")
    refs = [line.split(". ", 1)[1] for line in md.split("## References", 1)[1].splitlines()
            if line[:1].isdigit()]
    s += [PageBreak(), P("References", "h1")]
    for i, ref in enumerate(refs, 1):
        ref = ref.replace("&", "&amp;")
        ref = __import__("re").sub(r"\*([^*]+)\*", r"<i>\1</i>", ref)
        ref = __import__("re").sub(r"(https://doi\.org/\S+)", r"<link href='\1' color='#1f4e79'>\1</link>", ref)
        s.append(P(f"[{i}] {ref}", "ref"))
    s.append(P("<i>All references were checked against their Crossref records (authors, title, journal, volume, "
               "pages, year and DOI).</i>", "cap"))
    return s


def footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont(BODY, 8)
    canvas.setFillColor(colors.grey)
    if doc.page > 1:
        canvas.drawString(2 * cm, 1.2 * cm, "Convergent evolution of luciferases — final report")
        canvas.drawRightString(A4[0] - 2 * cm, 1.2 * cm, str(doc.page))
    canvas.restoreState()


def run() -> Path:
    data = load()
    make_figures(data)
    doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=1.8 * cm, bottomMargin=2 * cm,
                            title="Molecular Signatures of Convergent Evolution", author="Bioinformatics project")
    doc.build(build(data), onFirstPage=footer, onLaterPages=footer)
    print(f"wrote {OUT.relative_to(ROOT)}")
    return OUT


if __name__ == "__main__":
    run()

#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path

import markdown
from weasyprint import HTML


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "", "NA") else math.nan
    except ValueError:
        return math.nan


def fmt(value: str | float, digits: int = 3) -> str:
    x = num(str(value))
    if not math.isfinite(x):
        return "NA"
    if x != 0 and abs(x) < 0.001:
        return f"{x:.2e}"
    return f"{x:.{digits}f}"


def table(headers: list[str], rows: list[list[object]]) -> str:
    def clean(value: object) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")
    out = ["| " + " | ".join(map(clean, headers)) + " |",
           "| " + " | ".join(["---"] * len(headers)) + " |"]
    out += ["| " + " | ".join(clean(x) for x in row) + " |" for row in rows]
    return "\n".join(out)


def top(rows, predicate, n, key):
    out = [row for row in rows if predicate(row)]
    out.sort(key=key)
    return out[:n]


CSS = """
@page { size:A4; margin:16mm 15mm 18mm;
  @bottom-left { content:"GSE42568 reproducible analysis"; color:#64748b; font-size:8pt; }
  @bottom-right { content:"Page " counter(page) " of " counter(pages); color:#64748b; font-size:8pt; }}
@page:first { @bottom-left { content:none; } @bottom-right { content:none; }}
html { font-family:"DejaVu Sans",Arial,sans-serif; color:#0f172a; font-size:10pt; line-height:1.46; }
body { margin:0; }
.cover { min-height:245mm; display:flex; flex-direction:column; justify-content:center; page-break-after:always; padding:0 7mm; }
.eyebrow { color:#0e7490; font-weight:700; letter-spacing:.16em; text-transform:uppercase; }
.cover h1 { font-size:30pt; line-height:1.08; border:0; margin:7mm 0 5mm; }
.subtitle { color:#475569; font-size:14pt; max-width:155mm; }
.meta { border-top:2px solid #c2415d; color:#475569; margin-top:18mm; padding-top:5mm; }
h1 { font-size:23pt; border-bottom:3px solid #c2415d; padding-bottom:3mm; margin:0 0 7mm; }
h2 { font-size:16pt; margin:8mm 0 3mm; page-break-after:avoid; }
h3 { color:#0e7490; font-size:12pt; margin:5mm 0 2mm; page-break-after:avoid; }
p { margin:0 0 3.5mm; }
ul,ol { margin:2mm 0 4mm 6mm; padding-left:4mm; }
li { margin-bottom:1.2mm; }
a { color:#0e7490; text-decoration:none; }
blockquote { background:#ecfeff; border-left:4px solid #0e7490; color:#334155; margin:4mm 0; padding:4mm 5mm; }
code { background:#eef2f7; border-radius:1mm; font-size:8.5pt; padding:.4mm 1mm; }
pre { background:#0f172a; color:white; font-size:8pt; padding:4mm; white-space:pre-wrap; }
table { border-collapse:collapse; font-size:8pt; margin:3mm 0 6mm; width:100%; }
thead { display:table-header-group; } tr { page-break-inside:avoid; }
th { background:#e2e8f0; text-align:left; } th,td { border-bottom:.4pt solid #cbd5e1; padding:2mm; vertical-align:top; }
tbody tr:nth-child(even) { background:#f8fafc; }
img { display:block; margin:4mm auto 2mm; max-height:225mm; max-width:100%; object-fit:contain; }
p:has(> img) { page-break-inside:avoid; }
"""


def render(project: Path, stem: str, title: str, subtitle: str, label: str, body_md: str) -> None:
    reports = project / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (project / "output" / "pdf").mkdir(parents=True, exist_ok=True)
    (reports / f"{stem}.md").write_text(body_md.rstrip() + "\n", encoding="utf-8")
    cover = f'''<div class="cover"><div class="eyebrow">{html.escape(label)}</div><h1>{html.escape(title)}</h1>
<div class="subtitle">{html.escape(subtitle)}</div><div class="meta"><b>Dataset:</b> GEO GSE42568 · GPL570<br>
<b>Contrast:</b> cancer minus normal breast tissue<br><b>Release:</b> v1.0 · 11 September 2026<br>
<b>Scientific status:</b> PASS WITH DOCUMENTED LIMITATIONS</div></div>'''
    body = markdown.markdown(body_md, extensions=["tables", "fenced_code", "toc"])
    document = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>{CSS}</style></head><body>{cover}<main>{body}</main></body></html>'''
    html_path = reports / f"{stem}.html"
    html_path.write_text(document, encoding="utf-8")
    HTML(string=document, base_url=str(reports)).write_pdf(str(project / "output" / "pdf" / f"{stem}.pdf"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    project = parser.parse_args().project.resolve()
    tdir = project / "results" / "tables"
    de = read_tsv(tdir / "de_primary_all_genes.tsv")
    stability = read_tsv(tdir / "de_sensitivity_summary.tsv")
    loro = read_tsv(tdir / "leave_one_run_out_summary.tsv")
    treat = read_tsv(tdir / "treat_threshold_sensitivity_counts.tsv")
    go = read_tsv(tdir / "go_ora_all_terms.tsv")
    gsea = read_tsv(tdir / "gsea_all_gene_sets.tsv")
    annotations = read_tsv(tdir / "annotation_metrics.tsv")
    matrix = read_tsv(tdir / "matrix_integrity.tsv")
    review = read_tsv(tdir / "qc_candidate_review.tsv")
    licenses = read_tsv(tdir / "gene_set_manifest.tsv")
    acceptance = read_tsv(tdir / "acceptance_checks.tsv")

    discovered = sum(r["statistically_discovered"] == "TRUE" for r in de)
    priority = sum(r["priority_effect"] == "TRUE" for r in de)
    headline = sum(r["headline_eligible"] == "TRUE" for r in de)
    n_up = sum(r["priority_effect"] == "TRUE" and num(r["log2FC"]) > 0 for r in de)
    n_down = priority - n_up

    up = top(de, lambda r: r["headline_eligible"] == "TRUE" and num(r["log2FC"]) > 0, 10,
             lambda r: (num(r["treat_bh_fdr"]), -abs(num(r["log2FC"])), r["entrez_id"]))
    down = top(de, lambda r: r["headline_eligible"] == "TRUE" and num(r["log2FC"]) < 0, 10,
               lambda r: (num(r["treat_bh_fdr"]), -abs(num(r["log2FC"])), r["entrez_id"]))
    gene_rows = [[r["symbol_at_run"] or r["entrez_id"], r["entrez_id"], fmt(r["log2FC"], 2),
                  fmt(r["fold_change_cancer_over_normal"], 2), f"{fmt(r['ci95_low'], 2)} to {fmt(r['ci95_high'], 2)}",
                  fmt(r["treat_bh_fdr"], 2), r["representative_probe_id"],
                  f"{100 * num(r['constituent_same_direction_fraction']):.0f}%"] for r in up + down]

    s_names = ["P1_vs_S1_group_only", "P1_vs_S2_mixed_runs", "P1_vs_QC_candidate_exclusion",
               "P1_vs_S4_no_trend", "P1_vs_S4_group_aware", "P1_vs_median_collapse"]
    s_rows = []
    for name in s_names:
        r = next(x for x in stability if x["comparison"] == name)
        s_rows.append([name.replace("P1_vs_", ""), r["n_samples"], fmt(r["spearman_all_tested"]),
                       fmt(r["pearson_all_tested"]), f"{100 * num(r['p1_priority_sign_concordance']):.1f}%",
                       fmt(r["median_abs_effect_ratio_p1_priority"])])

    go_up = top(go, lambda r: r.get("ontology") == "BP" and r.get("direction") == "cancer_up" and num(r.get("bh_fdr")) < .05,
                6, lambda r: (num(r["bh_fdr"]), -num(r["hit_count"]), r["term_id"]))
    go_down = top(go, lambda r: r.get("ontology") == "BP" and r.get("direction") == "normal_up" and num(r.get("bh_fdr")) < .05,
                  6, lambda r: (num(r["bh_fdr"]), -num(r["hit_count"]), r["term_id"]))
    go_rows = [[r["direction"], r["term_name"], r["term_id"], r["hit_count"], r["intersected_set_size"], fmt(r["bh_fdr"], 2)]
               for r in go_up + go_down]
    rx = top(gsea, lambda r: r.get("collection") == "Reactome_preranked_GSEA_fallback" and num(r.get("bh_fdr")) < .05,
             12, lambda r: (num(r["bh_fdr"]), -abs(num(r["NES"])), r["term_id"]))
    rx_rows = [[r["direction"], r["term_name"], r["term_id"], fmt(r["NES"], 2), r["intersected_set_size"], fmt(r["bh_fdr"], 2)] for r in rx]

    client = f"""
[TOC]

# Executive summary

This analysis compares **104 breast cancer biopsies with 17 normal breast tissue samples** in GEO GSE42568. It uses the deposited log2 GCRMA matrix and a run-adjusted limma model. The target is an average **cancer-minus-normal bulk-tissue expression association**—not a causal effect, clinical classifier, prognosis model, or treatment-response model.

The workflow passed **32/32 Scientific QA checks** and 13/13 reproducibility guards. A clean rebuild from the checksum-locked official inputs reproduced **10/10 key result tables byte-for-byte**. Of **20,357 unique Entrez genes**, **{discovered:,}** reached moderated BH FDR < 0.05; **{priority:,}** met TREAT BH FDR < 0.05 while directly testing an absolute effect greater than 1.2-fold ({n_up:,} cancer-up; {n_down:,} normal-up). **{headline:,}** also passed probe and sensitivity interpretation checks.

> **Release decision:** suitable as an auditable portfolio/client-style discovery analysis with explicit limitations. It is not a diagnostic product or evidence that a pathway is activated, causal, or druggable.

# Data and boundary

The official [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) reports 104 pretreatment breast cancer biopsies and 17 normal breast tissues on GPL570 and cites Clarke et al. ([PMID 23740839](https://pubmed.ncbi.nlm.nih.gov/23740839/)). Normal-tissue donor/source status is not established, so this report never calls them healthy, adjacent, reduction-mammoplasty, or matched controls.

- Matrix SHA-256: `43b93f72f1d81838dcd2a5b6ccc3bc1875e46230dad0d7346236773a99e32b76`.
- 54,675 probe sets × 121 samples; finite deposited log2 GCRMA signals.
- No second log, re-normalization, imputation or ComBat was applied to the DE response.
- P1: `expression ~ processing_run_proxy + tissue_group`; normal reference.
- Thirteen run-proxy levels; six mixed levels containing 62 samples.

![Feature flow](../results/figures/F01_analysis_feature_flow.png)

# QC and design

Annotation retained 39,815 eligible probes and resolved 20,357 unique Entrez genes by an outcome-blind representative rule. Four samples triggered multi-category processed-data warnings. They show no gross matrix failure and have no raw-array corroboration, so all remain in P1; a separate 117-sample sensitivity measures influence.

![Primary PCA](../results/figures/F05_primary_PCA.png)

The design is full rank (14/14; residual df 107). Cancer-only run levels do not provide within-run group contrast once their fixed effects are included, explaining why P1 and the 62-sample mixed-run S2 estimates are identical.

![Design stability](../results/figures/F08_design_stability.png)

{table(["Comparison", "n", "Spearman", "Pearson", "Priority sign", "Median |effect| ratio"], s_rows)}

The design traffic light is **Green**. Removing the normal-heavy 2005-01-25 run retains Spearman 0.953 and every P1-priority direction. Excluding the four reviewed candidates retains every priority direction (Spearman 0.912).

# Differential expression

Positive log2FC means higher expression in cancer. The following genes are the fixed ten-smallest-TREAT-FDR set per direction, with stable tie-breaks; familiarity did not affect selection.

{table(["Gene", "Entrez", "log2FC", "Cancer/normal FC", "95% CI", "TREAT FDR", "Probe", "Probe concordance"], gene_rows)}

![Volcano plot](../results/figures/F10_volcano_plot.png)

![Priority heatmap](../results/figures/F11_priority_gene_heatmap.png)

Strong cancer-up signals include epithelial/adhesion-associated **DSP, EPCAM, ESRP1, SPINT2, and TFAP2A**. Prominent normal-up signals include **LYVE1, CIDEA, ACSM5, GPD1, ADH1C, and ACADL**. These are cohort-level bulk-tissue patterns and cannot distinguish tumor-cell regulation from tissue-composition shifts.

# Functional enrichment

GO ORA uses direction-specific priority genes and resource-specific tested-and-annotated backgrounds (GO BP: 15,588 genes). Reactome uses 10,178 tested genes with human Reactome annotation. GSEA begins with all 20,357 P1 moderated-t statistics.

{table(["Direction", "GO BP term", "ID", "Hits", "Set size", "BH FDR"], go_rows)}

![GO BP ORA](../results/figures/F13_GO_BP_ORA_dotplot.png)

Cancer-up ranks emphasize chromosome segregation, cell-cycle, chromatin and DNA-repair programs. Normal-up ranks emphasize fatty-acid metabolism, aerobic respiration and GPCR-related programs.

{table(["Direction", "Reactome set", "ID", "NES", "Set size", "BH FDR"], rx_rows)}

![Reactome GSEA](../results/figures/F14_Reactome_GSEA_NES.png)

![Selected GSEA curves](../results/figures/F15_selected_enrichment_curves.png)

Enrichment identifies concentration in the statistical ranking; it does not prove activation, causation or therapeutic vulnerability. Related GO terms overlap and are not independent findings.

# Licence decisions

- GO and Reactome completed with pinned sources and attribution; [Reactome data licence](https://reactome.org/license).
- KEGG: `KEGG_NOT_RUN_LICENSE_GATE`; no API, definition, cache or map used.
- MSigDB Hallmark: not downloaded or redistributed because commercial-use review was unconfirmed.
- Reactome is explicitly the licensed GSEA fallback.

# Limitations

1. Normal-tissue provenance, age and pairing identifiers are absent.
2. This is an unpaired, observational, cross-sectional bulk-tissue analysis; composition can drive signals.
3. Run is partially confounded with group even though the coefficient is estimable and sensitivities are Green.
4. Processed data cannot provide CEL-level NUSE/RLE, degradation, image-artifact or alternative-preprocessing checks.
5. Cancers span histology, grade and ER status; v1 intentionally estimates a broad average.
6. No external validation cohort was added to this bounded delivery.

# Deliverables and references

Complete local-workspace results are in `results/tables/`; F01–F15 each have PDF, PNG, SVG, caption and source TSV. The public portfolio edition intentionally omits large/generated source artifacts and has its own acceptance record in `docs/PUBLIC_PORTFOLIO_ACCEPTANCE.md`.

References: [GEO GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568); Clarke et al., *Carcinogenesis* 2013 ([doi:10.1093/carcin/bgt208](https://doi.org/10.1093/carcin/bgt208)); Ritchie et al., limma ([doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)); Subramanian et al., GSEA ([PMC1239896](https://pmc.ncbi.nlm.nih.gov/articles/PMC1239896/)); [GO citation policy](https://geneontology.org/docs/go-citation-policy/); [Reactome licence](https://reactome.org/license).
"""

    ann_rows = [[r["metric"], r["count"]] for r in annotations]
    treat_rows = [[r["minimum_fold_change"], fmt(r["minimum_abs_log2FC"]), r["treat_fdr_lt_0_05"], r["up"], r["down"]] for r in treat]
    license_rows = [[r["collection"], r["status"], r["source_version"], r["universe_size"], r["eligible_gene_sets"], r["ranking"]] for r in licenses]
    methods = f"""
[TOC]

# Scope and inputs

The formal contract is `GSE42568_PROJECT_DESIGN.md`: one cancer-minus-normal contrast with no WGCNA, survival, classifier, immune-infiltration or drug-prediction expansion. `scripts/download_inputs.sh` downloads the official matrix and GPL570 annotation snapshot. The matrix must match the frozen byte size, gzip test and SHA-256 before parsing.

The accepted matrix has 54,675 probe sets, 121 samples, range 2.3128–16.0469 and zero non-finite values. GEO identifies it as log2 GCRMA; no second log or normalization is performed.

# Metadata and design

The series-matrix header supplies original and normalized metadata in `sample_manifest.tsv`. With no pairing ID and no normal-sample ages, the design is unpaired and does not adjust age. Cancer-only clinical covariates cannot be included in a joint cancer-normal model.

For gene *g*, `E[y_gi] = beta_0g + processing_run_proxy_i + beta_group,g I(cancer_i)`. Normal is reference. The 121 × 14 design is rank 14, with 107 residual df, 13 run levels, six mixed levels and 62 mixed-run samples. `beta_group,g` is cancer-minus-normal log2FC.

# Processed-data QC

QC includes sample median/IQR/MAD/quantiles/range, Pearson and Spearman connectivity, primary PCA and robust PCA distance. Warnings use robust |z| >3. PCA uses top 5,000 pooled-MAD genes, centered but not scaled; top 1,000 and all-gene versions are sensitivity outputs. Robust distance uses up to 10 PCs reaching ≥50% variance and a 0.999 contour. A sample reaches review only with at least two independent warning categories; no automated exclusion is allowed.

# Annotation

{table(["Audit metric", "Count/value"], ann_rows)}

Mapping uses `hgu133plus2.db` 3.13.0 and `org.Hs.eg.db` 3.22.0. Stable keys are Entrez IDs; symbols are display snapshots. AFFX controls, constants, unmapped, one-to-many and `_x_at` probes are excluded. The representative order is `_at`, `_a_at`, `_s_at`; ties use pooled median then probe ID. Group labels and statistics are never used.

Median collapse and all eligible probe fits provide independent summarization/technical layers. A material opposite probe is opposite in direction and itself significant by probe-level 1.2-fold TREAT BH FDR.

# DE and sensitivity

P1 uses limma `lmFit`, `eBayes(trend=TRUE, robust=TRUE)`, cancer-minus-normal coefficient, BH across all 20,357 genes, and posterior-SE 95% CIs. Discovery is moderated BH FDR <0.05. Priority is `treat(lfc=log2(1.2))` BH FDR <0.05.

{table(["Minimum FC", "Minimum |log2FC|", "TREAT FDR <0.05", "Up", "Down"], treat_rows)}

- S1: group only, all 121 samples.
- S2: run + group in six mixed runs (n=62).
- S3: not applicable; no signed technical failures.
- QC influence: remove four reviewed candidates (n=117).
- S4a: `trend=FALSE`; S4b: group-aware `voomaByGroup` weights.
- Median collapse, probe-level and 13 leave-one-run-out fits complete the locked sensitivities.

{table(["Comparison", "n", "Spearman", "Pearson", "Priority sign", "Median |effect| ratio"], s_rows)}

Green requires Spearman ≥0.90 and priority-sign concordance ≥95% for both P1/S1 and P1/S2; Amber minima are 0.80. Final status is Green. Headline genes must retain direction across core models, agree with median collapse and have no material opposite probe.

# Visualization

F01–F15 are vector PDF/SVG and 600-dpi PNG with caption and source TSV in the complete analysis workspace. Adjusted PCA is visualization-only. Volcano categories use BH FDR despite a raw-P y-axis. Heatmap selection is TREAT-first, ≤25 genes per direction; row z-scores cap at ±2.5 and never enter DE.

# Enrichment and licensing

The universe is P1-tested unique Entrez genes annotated in the relevant resource. GO BP/CC/MF ORA tests cancer-up and normal-up separately, one-sided hypergeometric, post-intersection size 10–500, BH within ontology × direction. Reactome ORA uses the same rule. Reactome GSEA ranks all P1 genes by moderated t descending, stable Entrez tie-break, exponent 1, size 15–500 and `fgseaMultilevel`; observed exact ties are zero.

{table(["Collection", "Status", "Version", "Universe", "Eligible sets", "Ranking"], license_rows)}

# Runtime and reproducibility

The runtime is defined by the pinned Conda records in `config/`; commands resolve from the active `PATH` or an optional `BIOINFO_ENV_PREFIX`. The verified environment used R 4.5.3, limma 3.66.0, fgsea 1.36.2, hgu133plus2.db 3.13.0, org.Hs.eg.db 3.22.0, GO.db 3.22.0 and reactome.db 1.95.0. `tests/run_tests.R` contains 13 guards; `R/05_scientific_qa.R` contains 32 scientific acceptance checks. A clean rebuild reproduced 10/10 key machine-readable outputs with identical SHA-256 hashes.

Methods references: [limma](https://bioconductor.org/packages/limma), [hgu133plus2.db](https://bioconductor.org/packages/hgu133plus2.db), [GO](https://geneontology.org/docs/go-citation-policy/), [Reactome](https://reactome.org/license), [fgsea](https://bioconductor.org/packages/fgsea).
"""

    review_rows = [[r["gsm_id"], r["tissue_group"], r["processing_run_proxy"], r["independent_signal_count"],
                    r["distribution_flag"], r["correlation_flag"], r["pca_flag"]] for r in review]
    loro_rows = [[r["removed_run"], r["removed_normal_n"], r["removed_cancer_n"], fmt(r["spearman_all_tested"]),
                  f"{100 * num(r['p1_priority_sign_concordance']):.1f}%"] for r in loro]
    domain_rows = []
    for domain in sorted({r["domain"] for r in acceptance}):
        sub = [r for r in acceptance if r["domain"] == domain]
        domain_rows.append([domain, len(sub), sum(r["passed"] == "TRUE" for r in sub), ", ".join(r["status"] for r in sub)])
    qc = f"""
[TOC]

# Disposition

**PASS WITH DOCUMENTED LIMITATIONS:** 32/32 Scientific QA checks and 13/13 reproducibility guards pass; a clean rebuild reproduced 10/10 key tables byte-for-byte. All 121 samples enter P1. Four unusual samples remain after signed review and are removed only in an influence sensitivity.

# Input and metadata

{table(["Metric", "Observed"], [[r["metric"], r["value"]] for r in matrix])}

![Feature flow](../results/figures/F01_analysis_feature_flow.png)

The manifest has 104 cancer and 17 normal GPL570 samples. Normal ages, source category and pairing IDs are absent. GEO's narrative says 20 cancers under 50 while parsed sample fields contain 27 under 50 and 77 at least 50; the discrepancy is preserved, not overwritten.

![Metadata and run structure](../results/figures/F02_metadata_and_run_structure.png)

# Distribution and sample review

![Expression distributions](../results/figures/F03_expression_distributions.png)

![Sample QC metrics](../results/figures/F04_sample_qc_metrics.png)

{table(["Sample", "Group", "Run", "Signals", "Distribution", "Correlation", "PCA"], review_rows)}

The four profiles are unusual but finite and smooth, with plausible expression neighbours. No raw CEL metric corroborates technical failure. They are retained in P1 and jointly excluded only for influence assessment; formal S3 is not applicable.

# PCA and correlation

![Primary PCA](../results/figures/F05_primary_PCA.png)

![Visualization-only adjusted PCA](../results/figures/F06_visualization_only_adjusted_PCA.png)

![Sample correlation](../results/figures/F07_sample_correlation_heatmap.png)

PC1/PC2 explain 20.3%/10.9% of top-5,000-MAD variation. PC1 tissue R² is 0.678; tissue partial R² given run is 0.540 and run partial R² given tissue is 0.134. Adjusted PCA never enters inference.

# Annotation and statistical QC

{table(["Audit metric", "Count/value"], ann_rows)}

Eligibility is 73.03%, above the frozen 60% stop gate; 20,357 genes exceed the 10,000-gene gate. All 54,675 probes receive one exclusive status.

![Design stability](../results/figures/F08_design_stability.png)

{table(["Comparison", "n", "Spearman", "Pearson", "Priority sign", "Median |effect| ratio"], s_rows)}

{table(["Removed run", "Normal removed", "Cancer removed", "Spearman", "Priority sign"], loro_rows)}

The normal-heavy 2005-01-25 removal retains Spearman 0.953 and every priority direction. Median collapse is the least concordant major technical layer (rho 0.874; 97.0% priority direction), and its disagreements are filtered from headline claims rather than hidden.

# Enrichment and licence QA

GO BP universe is 15,588; Reactome universe 10,178. GSEA has all 20,357 moderated-t values exactly once and zero exact ties. KEGG and MSigDB Hallmark have explicit expected gate rows; Reactome is the completed fallback.

# Acceptance summary

{table(["Domain", "Checks", "Passed", "Statuses"], domain_rows)}

The full row-level ledger is `results/tables/acceptance_checks.tsv`. Residual limitations are normal provenance/age, unpaired bulk tissue, partial run confounding, cancer heterogeneity, processed-only QC and no external validation.
"""

    render(project, "client_report", "GSE42568 Breast Cancer Expression Analysis",
           "Run-adjusted differential expression, sensitivity analysis and licensed functional enrichment",
           "Client-style scientific report", client)
    render(project, "methods_appendix", "GSE42568 Methods Appendix",
           "Frozen parameters, statistical models, gene universes and reproducibility controls",
           "Technical specification and methods", methods)
    render(project, "qc_appendix", "GSE42568 QC Appendix",
           "Input, metadata, processed-expression, annotation, design and sensitivity review",
           "Scientific quality-control record", qc)

    readme = f"""# GSE42568 — breast cancer expression analysis

> **Portfolio case study**
>
> End-to-end analysis of a public breast-cancer expression dataset, from GEO data acquisition through QC, differential expression, pathway enrichment, reproducible reporting and scientific QA.

**Deliverables demonstrated:** QC · PCA · differential expression · volcano/heatmap · GO/Reactome enrichment · GSEA · reproducible code · publication-ready figures · client report

<table>
<tr>
<td width="50%"><img src="results/figures/F10_volcano_plot.png" alt="Differential-expression volcano plot"><br><sub>Differential-expression evidence</sub></td>
<td width="50%"><img src="results/figures/F11_priority_gene_heatmap.png" alt="Priority-gene expression heatmap"><br><sub>Priority-gene expression patterns</sub></td>
</tr>
<tr>
<td width="50%"><img src="results/figures/F13_GO_BP_ORA_dotplot.png" alt="GO Biological Process enrichment"><br><sub>GO Biological Process enrichment</sub></td>
<td width="50%"><img src="results/figures/F14_Reactome_GSEA_NES.png" alt="Reactome preranked GSEA"><br><sub>Reactome preranked GSEA</sub></td>
</tr>
</table>

## What this case demonstrates

This bounded client-style workflow analyzes **104 breast cancer** and **17 normal breast tissue** samples from [GEO GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568), GPL570. It turns an official processed expression matrix into an auditable delivery with sample and expression QC, current probe annotation, run-adjusted limma modelling, sensitivity analysis, licensed enrichment, client-facing reports and formal scientific QA.

The target is a cross-sectional **cancer-minus-normal bulk-tissue association**, not diagnosis, prognosis, causality or treatment response.

## Results snapshot

- 20,357 Entrez genes tested; {discovered:,} moderated BH discoveries.
- {priority:,} TREAT 1.2-fold priorities ({n_up:,} cancer-up; {n_down:,} normal-up); {headline:,} headline eligible.
- P1/S1 rho 0.936; P1/S2 rho 1.000; priority directions 100% concordant.
- Full analysis QA 32/32; reproducibility guards 13/13; clean-rebuild key hashes 10/10 identical; design status **Green**.
- KEGG and MSigDB content were not accessed; Reactome is the documented licensed fallback.

## Reproduce on Linux

The repository contains no machine-specific project or Conda path. From the repository root:

```bash
conda env create -f config/environment.yml
conda activate bioinfo
bash scripts/run_all.sh
```

The runner uses `Rscript`, `python` and `conda` from the active `PATH`. An explicit environment prefix remains available when desired:

```bash
BIOINFO_ENV_PREFIX="$CONDA_PREFIX" bash scripts/run_all.sh
```

Reproducibility is provided by the pinned Conda environment: `config/environment.yml` is the portable specification and `config/conda-explicit-linux-64.txt` is the exact verified Linux package record. Python freeze and R `sessionInfo()` are also included. No separate `renv.lock` is used.

The downloader enforces official GEO byte sizes, gzip integrity and SHA-256. Downloaded inputs and generated matrices are intentionally absent from the public edition and are reconstructed by the pipeline.

## Public-edition integrity

The complete analysis workspace passed 32/32 Scientific QA checks. The public portfolio edition intentionally omits data, run-specific logs, high-volume figure source tables and full-workspace inventories. Its separate inclusion and integrity contract is documented in [Public portfolio acceptance](docs/PUBLIC_PORTFOLIO_ACCEPTANCE.md); QA031 remains explicitly scoped to the complete analysis workspace.

## Navigate

- [Client report](reports/client_report.html) · [PDF](output/pdf/client_report.pdf)
- [Methods appendix](reports/methods_appendix.html) · [PDF](output/pdf/methods_appendix.pdf)
- [QC appendix](reports/qc_appendix.html) · [PDF](output/pdf/qc_appendix.pdf)
- [Scientific QA](docs/SCIENTIFIC_QA.md) · [Formal design](GSE42568_PROJECT_DESIGN.md)
- [Delivery acceptance and project map](docs/DELIVERY_ACCEPTANCE.md)
- [Priority DE results](results/tables/de_priority_genes.tsv) · [Gene-set/licence manifest](results/tables/gene_set_manifest.tsv)

## Limits

Normal-tissue provenance, age and pairing are unknown; run is partially confounded; bulk composition and cancer heterogeneity remain; processed data cannot provide CEL-level QC; v1 has no external validation. See [limitations](docs/limitations.md).

Code is MIT licensed. GEO/source attribution and GO/Reactome versions are retained; KEGG and MSigDB content is not redistributed.
"""
    (project / "README.md").write_text(readme, encoding="utf-8")
    (project / "docs" / "limitations.md").write_text("""# Limitations

1. Normal tissues are not sufficiently described as healthy, adjacent, reduction-mammoplasty, or paired controls.
2. Normal ages and pairing IDs are absent; age-adjusted and paired models are unsupported.
3. The cross-sectional bulk-tissue association is not causal; composition can drive results.
4. Run is partially confounded with group despite an estimable, Green-sensitivity contrast.
5. Processed data cannot support CEL-level NUSE/RLE, degradation, spatial artefact or alternative-CDF checks.
6. Cancers span histology, grade and ER status; v1 estimates a broad average.
7. Annotation and symbols reflect pinned package snapshots; Entrez is the stable key.
8. GO/Reactome sets overlap; enrichment does not prove activation or treatment response.
9. No external cohort or experimental validation is included.
10. KEGG and MSigDB Hallmark remain behind commercial-use licence review gates.
11. Family SOFT was not retained because the locked matrix header supplied required fields; the omission is recorded.
""", encoding="utf-8")
    (project / "docs" / "decision_log.md").write_text("""# Decision log

## 2026-09-11 — v1 locked implementation

- Accepted the checksum-locked GEO processed matrix; no second log or normalization.
- Locked P1 to run proxy plus tissue group, normal reference, cancer-minus-normal coefficient.
- Used current Entrez annotation and an outcome-blind representative-probe rule.
- Retained four reviewed QC candidates in P1; no raw evidence established failure. Removed them only in an influence model; formal S3 is not applicable.
- Completed group-only, mixed-run, no-trend, group-aware, median-collapse, probe-level and leave-one-run-out sensitivities; final traffic light Green.
- Kept TREAT 1.2-fold priority and 1.1/1.5 count sensitivities unchanged.
- Ran GO and Reactome with resource-specific tested-and-annotated universes.
- Did not access KEGG or MSigDB content because commercial-use review was not documented.
- Kept family SOFT as a documented infrastructure omission; scope was not expanded.
""", encoding="utf-8")
    print("[report] Wrote Markdown, HTML and PDF reports plus README and decision records")


if __name__ == "__main__":
    main()

# GSE42568 — breast cancer expression analysis

**English** ｜ [简体中文](README.zh-CN.md)

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

- 20,357 Entrez genes tested; 5,381 moderated BH discoveries.
- 2,835 TREAT 1.2-fold priorities (1,497 cancer-up; 1,338 normal-up); 2,745 headline eligible.
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

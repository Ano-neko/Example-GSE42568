
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

| Audit metric | Count/value |
| --- | --- |
| mapping_status:primary_eligible | 39815 |
| mapping_status:unmapped | 9953 |
| mapping_status:one_to_many | 1524 |
| mapping_status:constant | 96 |
| mapping_status:cross_hybridizing_x_at | 3225 |
| mapping_status:control_affx | 62 |
| suffix_class:s_at | 11292 |
| suffix_class:at | 37394 |
| suffix_class:a_at | 1720 |
| suffix_class:x_at | 4207 |
| suffix_class:control | 62 |
| total_probe_sets | 54675 |
| noncontrol_nonconstant_denominator | 54517 |
| primary_eligible_probes | 39815 |
| primary_eligible_fraction | 0.730322651649944 |
| unique_testable_entrez | 20357 |
| genes_with_multiple_eligible_probes | 10281 |

Mapping uses `hgu133plus2.db` 3.13.0 and `org.Hs.eg.db` 3.22.0. Stable keys are Entrez IDs; symbols are display snapshots. AFFX controls, constants, unmapped, one-to-many and `_x_at` probes are excluded. The representative order is `_at`, `_a_at`, `_s_at`; ties use pooled median then probe ID. Group labels and statistics are never used.

Median collapse and all eligible probe fits provide independent summarization/technical layers. A material opposite probe is opposite in direction and itself significant by probe-level 1.2-fold TREAT BH FDR.

# DE and sensitivity

P1 uses limma `lmFit`, `eBayes(trend=TRUE, robust=TRUE)`, cancer-minus-normal coefficient, BH across all 20,357 genes, and posterior-SE 95% CIs. Discovery is moderated BH FDR <0.05. Priority is `treat(lfc=log2(1.2))` BH FDR <0.05.

| Minimum FC | Minimum \|log2FC\| | TREAT FDR <0.05 | Up | Down |
| --- | --- | --- | --- | --- |
| 1.1 | 0.138 | 3730 | 2027 | 1703 |
| 1.2 | 0.263 | 2835 | 1497 | 1338 |
| 1.5 | 0.585 | 1451 | 680 | 771 |

- S1: group only, all 121 samples.
- S2: run + group in six mixed runs (n=62).
- S3: not applicable; no signed technical failures.
- QC influence: remove four reviewed candidates (n=117).
- S4a: `trend=FALSE`; S4b: group-aware `voomaByGroup` weights.
- Median collapse, probe-level and 13 leave-one-run-out fits complete the locked sensitivities.

| Comparison | n | Spearman | Pearson | Priority sign | Median \|effect\| ratio |
| --- | --- | --- | --- | --- | --- |
| S1_group_only | 121 | 0.936 | 0.975 | 100.0% | 1.069 |
| S2_mixed_runs | 62 | 1.000 | 1.000 | 100.0% | 1.000 |
| QC_candidate_exclusion | 117 | 0.912 | 0.975 | 100.0% | 1.100 |
| S4_no_trend | 121 | 1.000 | 1.000 | 100.0% | 1.000 |
| S4_group_aware | 121 | 0.962 | 0.993 | 100.0% | 1.006 |
| median_collapse | 121 | 0.874 | 0.902 | 97.0% | 0.960 |

Green requires Spearman ≥0.90 and priority-sign concordance ≥95% for both P1/S1 and P1/S2; Amber minima are 0.80. Final status is Green. Headline genes must retain direction across core models, agree with median collapse and have no material opposite probe.

# Visualization

F01–F15 are vector PDF/SVG and 600-dpi PNG with caption and source TSV in the complete analysis workspace. Adjusted PCA is visualization-only. Volcano categories use BH FDR despite a raw-P y-axis. Heatmap selection is TREAT-first, ≤25 genes per direction; row z-scores cap at ±2.5 and never enter DE.

# Enrichment and licensing

The universe is P1-tested unique Entrez genes annotated in the relevant resource. GO BP/CC/MF ORA tests cancer-up and normal-up separately, one-sided hypergeometric, post-intersection size 10–500, BH within ontology × direction. Reactome ORA uses the same rule. Reactome GSEA ranks all P1 genes by moderated t descending, stable Entrez tie-break, exponent 1, size 15–500 and `fgseaMultilevel`; observed exact ties are zero.

| Collection | Status | Version | Universe | Eligible sets | Ranking |
| --- | --- | --- | --- | --- | --- |
| GO_BP_ORA | completed | org.Hs.eg.db 3.22.0; GO.db 3.22.0 | 15588 | 6187 | NA |
| GO_CC_ORA | completed | org.Hs.eg.db 3.22.0; GO.db 3.22.0 | 16428 | 756 | NA |
| GO_MF_ORA | completed | org.Hs.eg.db 3.22.0; GO.db 3.22.0 | 16139 | 1233 | NA |
| Reactome_ORA | completed | 1.95.0 | 10178 | 1658 | NA |
| Reactome_preranked_GSEA_fallback | completed | 1.95.0 | 10178 | 1367 | P1 moderated t, descending; exact ties by Entrez ID |
| KEGG_ORA | not_run_license_gate | NA | NA | NA | NA |
| MSigDB_Hallmark_GSEA | not_run_license_review_unconfirmed | NA | NA | NA | NA |

# Runtime and reproducibility

The runtime is defined by the pinned Conda records in `config/`; commands resolve from the active `PATH` or an optional `BIOINFO_ENV_PREFIX`. The verified environment used R 4.5.3, limma 3.66.0, fgsea 1.36.2, hgu133plus2.db 3.13.0, org.Hs.eg.db 3.22.0, GO.db 3.22.0 and reactome.db 1.95.0. `tests/run_tests.R` contains 13 guards; `R/05_scientific_qa.R` contains 32 scientific acceptance checks. A clean rebuild reproduced 10/10 key machine-readable outputs with identical SHA-256 hashes.

Methods references: [limma](https://bioconductor.org/packages/limma), [hgu133plus2.db](https://bioconductor.org/packages/hgu133plus2.db), [GO](https://geneontology.org/docs/go-citation-policy/), [Reactome](https://reactome.org/license), [fgsea](https://bioconductor.org/packages/fgsea).

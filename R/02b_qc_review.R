source("R/common.R")
ctx <- init_project()
ensure_dirs()

message("[02b] Applying signed Scientific Architect QC review")
ledger_path <- "results/tables/sample_exclusion_ledger.tsv"
metrics_path <- "results/tables/sample_qc_metrics.tsv"
ledger <- fread(ledger_path)
metrics <- fread(metrics_path)

expected <- c("GSM1045192", "GSM1045193", "GSM1045217", "GSM1045302")
observed <- sort(metrics[candidate_exclusion_review == TRUE, gsm_id])
stopifnot(identical(sort(expected), observed))
stopifnot(identical(sort(ledger$gsm_id), observed))
stopifnot(all(ledger$raw_evidence_flag == FALSE))

ledger[, `:=`(
  decision = "retain_primary_exclude_in_qc_candidate_sensitivity",
  rationale = paste(
    "Retained in the primary all-sample analysis because deposited-expression distributions are continuous and biologically plausible,",
    "no gross array failure is visible, and no raw-array evidence is available. Excluded jointly only in the locked QC-candidate sensitivity analysis."
  ),
  reviewer = "Scientific Architect",
  review_date = "2026-09-11",
  primary_analysis_included = TRUE,
  qc_candidate_sensitivity_included = FALSE
)]

evidence_cols <- c(
  "gsm_id", "tissue_group", "processing_run_proxy", "sample_median", "sample_iqr", "sample_mad",
  "median_pearson_correlation", "robust_pca_distance", "robust_pca_cutoff",
  "distribution_flag", "correlation_flag", "pca_flag", "raw_evidence_flag", "independent_signal_count"
)
evidence <- metrics[gsm_id %in% expected, ..evidence_cols]
evidence <- merge(evidence, ledger[, .(gsm_id, decision, rationale, reviewer, review_date,
                                      primary_analysis_included, qc_candidate_sensitivity_included)],
                  by = "gsm_id", all.x = TRUE, sort = FALSE)

write_tsv(ledger, ledger_path)
write_tsv(evidence, "results/tables/qc_candidate_review.tsv")

doc <- c(
  "# QC candidate review",
  "",
  "**Status:** reviewed and released for primary analysis  ",
  "**Reviewer:** Scientific Architect  ",
  "**Review date:** 2026-09-11",
  "",
  "## Trigger",
  "",
  "Four samples triggered at least two independent processed-expression warning categories. This was a pre-specified pause/review condition, not an automatic exclusion rule.",
  "",
  "## Evidence considered",
  "",
  "- All four arrays have smooth, finite expression distributions within the continuous range of the study; none shows a grossly truncated, singular, or globally corrupted profile.",
  "- Pearson connectivity and robust PCA confirm that the samples are unusual, especially GSM1045193, GSM1045217, and GSM1045302, but each retains plausible biological nearest neighbours.",
  "- GSM1045192 is supported by distribution and PCA warnings but does not cross the correlation warning boundary.",
  "- The deposited matrix provides no raw CEL quality metrics (for example RNA degradation, NUSE/RLE, or image artefacts) that could independently establish technical failure.",
  "- The warning pattern is not sufficient to distinguish technical quality from true tissue heterogeneity, particularly because cancer and normal breast tissue differ substantially and run is partially confounded with tissue group.",
  "",
  "## Decision",
  "",
  "All 121 samples are retained in the primary model. GSM1045192, GSM1045193, GSM1045217, and GSM1045302 are excluded jointly in the locked QC-candidate sensitivity model. This deliberately tests influence without retroactively labelling any sample as a failed array.",
  "",
  "A primary conclusion is considered robust only if its direction and magnitude remain materially consistent in this and the other pre-specified sensitivity analyses. The machine-readable signed decision is in `results/tables/sample_exclusion_ledger.tsv`; the supporting values are in `results/tables/qc_candidate_review.tsv`.",
  ""
)
writeLines(doc, "docs/qc_review.md")

message("[02b] PASS: all 121 samples retained in primary; 4 candidates locked for sensitivity exclusion")

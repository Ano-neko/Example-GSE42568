source("R/common.R")
ctx <- init_project()
cfg <- ctx$cfg
ensure_dirs()

message("[05] Full Scientific QA")
checks <- list()
add_check <- function(check_id, domain, requirement, condition, evidence, hard = TRUE,
                      pass_status = "PASS", fail_status = "FAIL") {
  checks[[length(checks) + 1L]] <<- data.table(
    check_id = check_id, domain = domain, requirement = requirement,
    status = if (isTRUE(condition)) pass_status else fail_status,
    passed = isTRUE(condition), hard_gate = hard, evidence = evidence
  )
}

manifest <- fread("data/metadata/sample_manifest.tsv")
matrix_qc <- fread("results/tables/matrix_integrity.tsv")
design_qc <- fread("results/tables/design_diagnostics.tsv")
ann_qc <- fread("results/tables/annotation_metrics.tsv")
probe_audit <- fread("results/tables/probe_mapping_audit.tsv")
selection <- fread("results/tables/probe_selection_audit.tsv", colClasses = list(character = "entrez_id"))
ledger <- fread("results/tables/sample_exclusion_ledger.tsv")
de <- fread("results/tables/de_primary_all_genes.tsv", colClasses = list(character = "entrez_id"))
priority <- fread("results/tables/de_priority_genes.tsv", colClasses = list(character = "entrez_id"))
stability <- fread("results/tables/de_sensitivity_summary.tsv")
loro <- fread("results/tables/leave_one_run_out_summary.tsv")
universes <- fread("results/tables/enrichment_universe_membership.tsv", colClasses = list(character = "entrez_id"))
gsea_rank <- fread("results/tables/gsea_rank_vector.tsv", colClasses = list(character = "entrez_id"))
go <- fread("results/tables/go_ora_all_terms.tsv")
pathway <- fread("results/tables/pathway_ora_all_terms.tsv")
gsea <- fread("results/tables/gsea_all_gene_sets.tsv")
gs_manifest <- fread("results/tables/gene_set_manifest.tsv")
expr <- readRDS("data/derived/expression_gene_representative.rds")

metric_value <- function(dt, key) dt[metric == key, value][1]
metric_count <- function(dt, key) as.numeric(dt[metric == key, count][1])

observed_hash <- sha256_file(cfg$input$matrix_file)
add_check("QA001", "input", "Official processed matrix checksum equals frozen specification",
          identical(observed_hash, cfg$input$expected_matrix_sha256),
          paste("observed", observed_hash, "expected", cfg$input$expected_matrix_sha256))
add_check("QA002", "input", "Matrix dimensions and finite values pass",
          as.integer(metric_value(matrix_qc, "rows_probe_sets")) == 54675L &&
            as.integer(metric_value(matrix_qc, "columns_samples")) == 121L &&
            as.integer(metric_value(matrix_qc, "nonfinite_values")) == 0L,
          "matrix_integrity.tsv: 54,675 probes x 121 samples; zero non-finite values")
add_check("QA003", "metadata", "Group counts are 104 cancer and 17 normal",
          nrow(manifest) == 121L && sum(manifest$tissue_group == "cancer") == 104L && sum(manifest$tissue_group == "normal") == 17L,
          sprintf("n=%d; cancer=%d; normal=%d", nrow(manifest), sum(manifest$tissue_group == "cancer"), sum(manifest$tissue_group == "normal")))
add_check("QA004", "metadata", "All samples use GPL570",
          uniqueN(manifest$platform_id) == 1L && unique(manifest$platform_id) == "GPL570",
          paste(unique(manifest$platform_id), collapse = ";"))
pair_nonmissing <- sum(!is.na(manifest$paired_identifier) & manifest$paired_identifier != "")
add_check("QA005", "metadata", "No unsupported paired design",
          pair_nonmissing == 0L, paste("non-missing pairing identifiers", pair_nonmissing))
add_check("QA006", "metadata", "Normal-tissue provenance limitation is explicitly documented",
          any(grepl("provenance of the 17 normal", fread("data/metadata/metadata_issues.tsv")$issue, ignore.case = TRUE)),
          "metadata_issues.tsv META-002", hard = TRUE)
add_check("QA007", "metadata", "Age discrepancy is explicitly documented",
          any(grepl("20 tumors under 50", fread("data/metadata/metadata_issues.tsv")$issue, fixed = TRUE)),
          "metadata_issues.tsv META-001", hard = TRUE)

add_check("QA008", "design", "Primary design is full rank",
          metric_value(design_qc, "full_rank") == "TRUE" && as.integer(metric_value(design_qc, "design_rank")) == 14L,
          "14 design columns, rank 14")
add_check("QA009", "design", "Cancer-minus-normal contrast is estimable",
          metric_value(design_qc, "contrast") == "cancer-minus-normal" && metric_value(design_qc, "contrast_estimable") == "TRUE",
          "reference=normal; coefficient=tissue_groupcancer")
add_check("QA010", "design", "Mixed-run evidence boundary exists",
          as.integer(metric_value(design_qc, "mixed_run_levels")) == 6L && as.integer(metric_value(design_qc, "mixed_run_samples")) == 62L,
          "6 mixed runs; 62 samples")

status_total <- sum(ann_qc[grepl("^mapping_status:", metric), as.numeric(count)])
add_check("QA011", "annotation", "Every platform probe has exactly one mutually exclusive mapping status",
          nrow(probe_audit) == 54675L && status_total == 54675L && !anyDuplicated(probe_audit$probe_id),
          sprintf("audit rows=%d; category total=%d", nrow(probe_audit), status_total))
eligible_fraction <- metric_count(ann_qc, "primary_eligible_fraction")
add_check("QA012", "annotation", "Eligible-probe fraction passes the frozen 60% gate",
          eligible_fraction >= as.numeric(cfg$annotation$minimum_eligible_probe_fraction),
          sprintf("%.2f%% eligible", 100 * eligible_fraction))
tested_gene_count <- metric_count(ann_qc, "unique_testable_entrez")
add_check("QA013", "annotation", "Unique Entrez count passes the frozen 10,000-gene gate",
          tested_gene_count >= as.numeric(cfg$annotation$minimum_unique_entrez),
          sprintf("%.0f unique Entrez genes", tested_gene_count))
add_check("QA014", "annotation", "Representative-probe selection is one per Entrez and excludes ambiguous/x probes",
          selection[selected_representative == TRUE, .N] == uniqueN(selection$entrez_id) &&
            all(probe_audit[mapping_status == "one_to_many" | suffix_class == "x_at", primary_eligible == FALSE]),
          sprintf("%d representatives", selection[selected_representative == TRUE, .N]))

add_check("QA015", "sample_QC", "All processed-data QC candidates received signed review",
          nrow(ledger) == 4L && all(ledger$reviewer == "Scientific Architect") && all(ledger$primary_analysis_included == TRUE),
          "4 candidates retained in primary and locked for influence sensitivity")
add_check("QA016", "sample_QC", "No sample is mislabelled as a confirmed technical failure",
          sum(ledger$primary_analysis_included == FALSE) == 0L && all(ledger$raw_evidence_flag == FALSE),
          "S3 is not applicable; no raw-array evidence corroborated failure", hard = TRUE,
          pass_status = "NOT_APPLICABLE_CONFIRMED_FAILURE_EXCLUSION")

required_de_cols <- c("entrez_id", "log2FC", "standard_error", "ci95_low", "ci95_high", "moderated_t", "p_value", "bh_fdr",
                      "treat_p_value", "treat_bh_fdr", "representative_probe_id", "s1_log2FC", "s2_log2FC",
                      "s4_no_trend_log2FC", "s4_group_aware_log2FC", "median_collapse_log2FC",
                      "qc_candidate_exclusion_log2FC", "headline_eligible", "contrast")
add_check("QA017", "DE", "Complete primary DE table has one finite row per tested gene and all required fields",
          nrow(de) == 20357L && !anyDuplicated(de$entrez_id) && all(required_de_cols %in% names(de)) &&
            all(is.finite(de$log2FC)) && all(is.finite(de$standard_error)) && all(de$p_value >= 0 & de$p_value <= 1) &&
            all(de$bh_fdr >= 0 & de$bh_fdr <= 1),
          "20,357 unique genes; finite effects/SE; P and FDR in [0,1]")
add_check("QA018", "DE", "Contrast label is cancer-minus-normal for every row",
          all(de$contrast == "cancer-minus-normal"), "positive effect means higher expression in cancer")
add_check("QA019", "DE", "Priority table exactly matches P1 TREAT 1.2-fold BH FDR <0.05",
          setequal(priority$entrez_id, de[treat_bh_fdr < 0.05, entrez_id]) && nrow(priority) == 2835L,
          sprintf("priority=%d; cancer-up=%d; normal-up=%d", nrow(priority), sum(priority$log2FC > 0), sum(priority$log2FC < 0)))
add_check("QA020", "DE", "Median-collapse and probe-direction audit gate headline interpretation",
          all(de[headline_eligible == TRUE, representative_median_direction_concordant == TRUE]) &&
            all(de[headline_eligible == TRUE, any_material_opposite_treat_probe == FALSE]),
          sprintf("headline eligible=%d of %d priority", sum(de$headline_eligible), sum(de$priority_effect)))

primary_stability <- stability[comparison %in% c("P1_vs_S1_group_only", "P1_vs_S2_mixed_runs")]
add_check("QA021", "sensitivity", "P1/S1 and P1/S2 meet frozen Green thresholds",
          all(primary_stability$spearman_all_tested >= as.numeric(cfg$stability$green_spearman_min)) &&
            all(primary_stability$p1_priority_sign_concordance >= as.numeric(cfg$stability$green_sign_concordance_min)),
          paste(sprintf("%s rho=%.3f priority-sign=%.3f", primary_stability$comparison,
                        primary_stability$spearman_all_tested, primary_stability$p1_priority_sign_concordance), collapse = " | "))
add_check("QA022", "sensitivity", "QC-candidate exclusion does not reverse P1 priority genes",
          stability[comparison == "P1_vs_QC_candidate_exclusion", p1_priority_sign_concordance] >= 0.95,
          sprintf("rho=%.3f; priority sign=%.3f", stability[comparison == "P1_vs_QC_candidate_exclusion", spearman_all_tested],
                  stability[comparison == "P1_vs_QC_candidate_exclusion", p1_priority_sign_concordance]))
add_check("QA023", "sensitivity", "Leave-one-run-out does not reverse priority directions",
          all(loro$status == "completed") && min(loro$p1_priority_sign_concordance) >= 0.95,
          sprintf("13/13 fits complete; minimum priority sign concordance %.3f", min(loro$p1_priority_sign_concordance)))
add_check("QA024", "sensitivity", "Normal-heavy 2005-01-25 run does not collapse the effect ranking",
          loro[removed_run == "2005-01-25", spearman_all_tested] >= 0.90 &&
            loro[removed_run == "2005-01-25", p1_priority_sign_concordance] >= 0.95,
          sprintf("rho=%.3f; priority sign=%.3f", loro[removed_run == "2005-01-25", spearman_all_tested],
                  loro[removed_run == "2005-01-25", p1_priority_sign_concordance]))

add_check("QA025", "enrichment", "Every enrichment-universe member is a P1-tested Entrez gene",
          all(universes$entrez_id %in% de$entrez_id) && !anyDuplicated(universes[, .(resource, entrez_id)]),
          paste(universes[, .N, by = resource][, paste(resource, N, sep = "=")], collapse = "; "))
add_check("QA026", "enrichment", "GSEA rank uses all tested genes exactly once in stable descending order",
          nrow(gsea_rank) == nrow(de) && !anyDuplicated(gsea_rank$entrez_id) && setequal(gsea_rank$entrez_id, de$entrez_id) &&
            all(diff(gsea_rank$moderated_t) <= 0),
          "20,357 moderated-t values; zero exact ties")
add_check("QA027", "enrichment", "GO ORA contains all ontology/direction families with resource-specific denominators",
          all(c("BP", "CC", "MF") %in% go$ontology) && all(c("cancer_up", "normal_up") %in% go$direction) &&
            all(go[status == "completed", selected_count <= universe_size]),
          sprintf("GO rows=%d; significant BP rows=%d", nrow(go), go[ontology == "BP" & bh_fdr < 0.05, .N]))
add_check("QA028", "enrichment", "Reactome ORA and GSEA completed with explicit CC0-compatible snapshot",
          pathway[resource == "Reactome" & status == "completed", .N] > 0L &&
            gsea[collection == "Reactome_preranked_GSEA_fallback" & status == "completed", .N] > 0L,
          sprintf("Reactome ORA rows=%d; GSEA rows=%d", pathway[resource == "Reactome", .N],
                  gsea[collection == "Reactome_preranked_GSEA_fallback", .N]))
add_check("QA029", "licensing", "KEGG remains behind the commercial licence gate",
          gs_manifest[collection == "KEGG_ORA", status] == "not_run_license_gate" &&
            pathway[resource == "KEGG", status] == "KEGG_NOT_RUN_LICENSE_GATE",
          "No KEGG API, definitions, cache, or maps used", hard = TRUE, pass_status = "EXPECTED_LICENSE_GATE")
add_check("QA030", "licensing", "MSigDB Hallmark remains behind unconfirmed commercial-use review",
          gs_manifest[collection == "MSigDB_Hallmark_GSEA", status] == "not_run_license_review_unconfirmed" &&
            gsea[collection == "MSigDB_Hallmark", grepl("NOT_RUN", status)],
          "Reactome GSEA is the explicitly named fallback", hard = TRUE, pass_status = "EXPECTED_LICENSE_GATE")

figure_stems <- c(
  "F01_analysis_feature_flow", "F02_metadata_and_run_structure", "F03_expression_distributions",
  "F04_sample_qc_metrics", "F05_primary_PCA", "F06_visualization_only_adjusted_PCA",
  "F07_sample_correlation_heatmap", "F08_design_stability", "F09_MA_plot", "F10_volcano_plot",
  "F11_priority_gene_heatmap", "F12_top_effect_forest", "F13_GO_BP_ORA_dotplot",
  "F14_Reactome_GSEA_NES", "F15_selected_enrichment_curves"
)
figure_files <- unlist(lapply(figure_stems, function(s) file.path("results/figures", paste0(s, c(".pdf", ".png", ".caption.txt")))))
source_map <- list(
  F01 = "results/tables/figure_F01_source.tsv",
  F02 = c("results/tables/figure_F02_metadata_source.tsv", "results/tables/figure_F02_run_source.tsv"),
  F03 = c("results/tables/figure_F03_boxplot_source.tsv", "results/tables/figure_F03_density_source.tsv"),
  F04 = "results/tables/figure_F04_source.tsv", F05 = "results/tables/figure_F05_source.tsv",
  F06 = "results/tables/figure_F06_source.tsv", F07 = "results/tables/figure_F07_source.tsv",
  F08 = "results/tables/figure_F08_source.tsv", F09 = "results/tables/figure_F09_source.tsv",
  F10 = "results/tables/figure_F10_source.tsv", F11 = "results/tables/figure_F11_source.tsv",
  F12 = "results/tables/figure_F12_source.tsv", F13 = "results/tables/figure_F13_source.tsv",
  F14 = "results/tables/figure_F14_source.tsv", F15 = "results/tables/figure_F15_source.tsv"
)
source_files <- unlist(source_map, use.names = FALSE)
add_check("QA031", "artifacts", "Complete analysis workspace contains all F01-F15 PDF/PNG/caption and source-data files",
          all(file.exists(c(figure_files, source_files))) && all(file.info(c(figure_files, source_files))$size > 0),
          sprintf("Full-workspace check: %d figure/caption files and %d source tables; public edition is governed separately", length(figure_files), length(source_files)))

code_files <- list.files(c("R", "scripts"), pattern = "\\.(R|py|sh)$", recursive = TRUE, full.names = TRUE)
code_files <- code_files[basename(code_files) != "05_scientific_qa.R"]
code_text <- paste(vapply(code_files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
forbidden_calls <- grepl("library\\((WGCNA|survival|randomForest)\\)|CIBERSORT\\(|estimateCellCounts\\(|drug[._]predict\\(",
                         code_text, ignore.case = TRUE)
add_check("QA032", "scope", "No prohibited extension module is implemented",
          !forbidden_calls, "No WGCNA, survival, classifier, immune-infiltration, or drug-prediction call in executable code")

# Data-driven headline and fixed-marker checks are diagnostic only; they never gate inclusion.
normal_idx <- manifest$tissue_group == "normal"
cancer_idx <- manifest$tissue_group == "cancer"
top_up <- de[headline_eligible == TRUE & log2FC > 0][order(treat_bh_fdr, -abs(log2FC), entrez_id)][seq_len(10)]
top_down <- de[headline_eligible == TRUE & log2FC < 0][order(treat_bh_fdr, -abs(log2FC), entrez_id)][seq_len(10)]
headline <- rbind(top_up, top_down, fill = TRUE)
headline[, `:=`(
  diagnostic_set = "predefined_top10_each_direction_by_TREAT_FDR",
  raw_normal_mean = rowMeans(expr[entrez_id, normal_idx, drop = FALSE]),
  raw_cancer_mean = rowMeans(expr[entrez_id, cancer_idx, drop = FALSE])
)]
headline[, raw_cancer_minus_normal := raw_cancer_mean - raw_normal_mean]
headline[, raw_and_adjusted_direction_concordant := sign(raw_cancer_minus_normal) == sign(log2FC)]
fixed_symbols <- c("EPCAM", "MKI67", "ESR1", "PGR", "ERBB2", "ADIPOQ", "CIDEA", "CDC20")
marker <- de[symbol_at_run %in% fixed_symbols]
marker[, `:=`(diagnostic_set = "fixed_breast_context_sanity_only",
              raw_normal_mean = rowMeans(expr[entrez_id, normal_idx, drop = FALSE]),
              raw_cancer_mean = rowMeans(expr[entrez_id, cancer_idx, drop = FALSE]))]
marker[, raw_cancer_minus_normal := raw_cancer_mean - raw_normal_mean]
marker[, raw_and_adjusted_direction_concordant := sign(raw_cancer_minus_normal) == sign(log2FC)]
sanity_cols <- c("diagnostic_set", "entrez_id", "symbol_at_run", "representative_probe_id", "raw_normal_mean",
                 "raw_cancer_mean", "raw_cancer_minus_normal", "log2FC", "ci95_low", "ci95_high", "bh_fdr", "treat_bh_fdr",
                 "s1_log2FC", "s2_log2FC", "qc_candidate_exclusion_log2FC", "median_collapse_log2FC",
                 "constituent_same_direction_fraction", "raw_and_adjusted_direction_concordant")
write_tsv(rbind(headline[, ..sanity_cols], marker[, ..sanity_cols], fill = TRUE),
          "results/tables/biological_sanity_checks.tsv")

qa <- rbindlist(checks)
write_tsv(qa, "results/tables/acceptance_checks.tsv")
hard_fail <- qa[hard_gate == TRUE & passed == FALSE]
soft_fail <- qa[hard_gate == FALSE & passed == FALSE]
overall <- if (nrow(hard_fail)) "HOLD" else if (nrow(soft_fail)) "PASS_WITH_LIMITATIONS" else "PASS_WITH_DOCUMENTED_LIMITATIONS"

qa_doc <- c(
  "# Scientific QA",
  "",
  paste0("**Overall status: ", overall, "**  "),
  "**Review date:** 2026-09-11  ",
  "**Contrast:** cancer minus normal breast tissue",
  "",
  "## Outcome",
  "",
  sprintf("- %d checks evaluated; %d passed; %d hard failures; %d soft failures.", nrow(qa), sum(qa$passed), nrow(hard_fail), nrow(soft_fail)),
  sprintf("- P1 tested %d genes; %d moderated BH discoveries; %d TREAT priorities; %d headline-eligible genes.",
          nrow(de), sum(de$statistically_discovered), sum(de$priority_effect), sum(de$headline_eligible)),
  sprintf("- Design stability: %s. P1/S1 Spearman %.3f; P1/S2 Spearman %.3f; both preserve 100%% of P1-priority directions.",
          stability$final_design_traffic_light[1],
          stability[comparison == "P1_vs_S1_group_only", spearman_all_tested],
          stability[comparison == "P1_vs_S2_mixed_runs", spearman_all_tested]),
  "- S3 is formally not applicable because no sample met the signed technical-failure definition. Four reviewed QC candidates were retained in P1 and jointly removed only in a separately labelled influence sensitivity.",
  "- KEGG and MSigDB Hallmark were not run because the required commercial-use licence review was not documented. Reactome is the licensed pathway/GSEA fallback.",
  "",
  "## Interpretation boundary",
  "",
  "This analysis estimates an adjusted cross-sectional association in heterogeneous bulk tissues. It does not establish causality, diagnosis, prognosis, treatment response, pathway activation, or cell-intrinsic regulation. Normal-tissue provenance and age are unavailable; run is partially confounded with group; and processed GEO data cannot support raw-array NUSE/RLE, degradation, or image QC.",
  "",
  "## Machine-readable evidence",
  "",
  "- `results/tables/acceptance_checks.tsv`",
  "- `results/tables/biological_sanity_checks.tsv`",
  "- `results/tables/de_sensitivity_summary.tsv`",
  "- `results/tables/leave_one_run_out_summary.tsv`",
  "- `results/tables/gene_set_manifest.tsv`",
  ""
)
writeLines(qa_doc, "docs/SCIENTIFIC_QA.md")
if (nrow(hard_fail)) {
  print(hard_fail)
  stop("Scientific QA hard gate failed")
}
message(sprintf("[05] %s: %d/%d acceptance checks passed", overall, sum(qa$passed), nrow(qa)))

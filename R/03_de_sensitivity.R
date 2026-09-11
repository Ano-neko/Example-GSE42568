source("R/common.R")
ctx <- init_project()
cfg <- ctx$cfg
ensure_dirs()

suppressPackageStartupMessages({
  library(pheatmap)
  library(parallel)
})

message("[03] Locked differential-expression and sensitivity analysis")
expr_rep <- readRDS("data/derived/expression_gene_representative.rds")
expr_med <- readRDS("data/derived/expression_gene_median_collapse.rds")
expr_probe_all <- readRDS("data/derived/expression_probe_matrix.rds")
pairs <- as.data.table(readRDS("data/derived/eligible_probe_pairs.rds"))
gene_ann <- as.data.table(readRDS("data/derived/gene_annotation.rds"))
metadata <- as.data.table(readRDS("data/derived/sample_metadata.rds"))
ledger <- fread("results/tables/sample_exclusion_ledger.tsv")

stopifnot(nrow(expr_rep) == nrow(gene_ann), nrow(expr_rep) == 20357L)
stopifnot(identical(rownames(expr_rep), rownames(expr_med)))
stopifnot(identical(colnames(expr_rep), metadata$gsm_id))
stopifnot(!anyDuplicated(rownames(expr_rep)), all(is.finite(expr_rep)), all(is.finite(expr_med)))
stopifnot(identical(sort(pairs$probe_id), sort(intersect(pairs$probe_id, rownames(expr_probe_all)))))

bh_cut <- as.numeric(cfg$differential_expression$bh_fdr)
fc_main <- as.numeric(cfg$differential_expression$treat_fold_change)
fc_sens <- as.numeric(unlist(cfg$differential_expression$treat_sensitivity_fold_changes))
lfc_main <- log2(fc_main)

build_design <- function(md, include_run = TRUE) {
  md <- as.data.frame(md)
  md$tissue_group <- factor(md$tissue_group, levels = c("normal", "cancer"))
  if (include_run) {
    md$processing_run_proxy <- factor(md$processing_run_proxy)
    design <- model.matrix(~ processing_run_proxy + tissue_group, data = md)
  } else {
    design <- model.matrix(~ tissue_group, data = md)
  }
  stopifnot(qr(design)$rank == ncol(design), "tissue_groupcancer" %in% colnames(design))
  design
}

extract_extended <- function(fit0, fit, treat_fit, coef_index) {
  se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = fit$df.total)
  data.table(
    feature_id = rownames(fit$coefficients),
    log2FC = fit$coefficients[, coef_index],
    fold_change_cancer_over_normal = 2 ^ fit$coefficients[, coef_index],
    average_expression = fit$Amean,
    standard_error = se,
    ci95_low = fit$coefficients[, coef_index] - crit * se,
    ci95_high = fit$coefficients[, coef_index] + crit * se,
    moderated_t = fit$t[, coef_index],
    p_value = fit$p.value[, coef_index],
    bh_fdr = p.adjust(fit$p.value[, coef_index], method = "BH"),
    treat_t = treat_fit$t[, coef_index],
    treat_p_value = treat_fit$p.value[, coef_index],
    treat_bh_fdr = p.adjust(treat_fit$p.value[, coef_index], method = "BH")
  )
}

fit_standard <- function(mat, md, include_run = TRUE, trend = TRUE, robust = TRUE,
                         treat_fc = fc_main) {
  design <- build_design(md, include_run)
  j <- match("tissue_groupcancer", colnames(design))
  fit0 <- lmFit(mat, design)
  fit <- eBayes(fit0, trend = trend, robust = robust)
  tr <- treat(fit0, lfc = log2(treat_fc), trend = trend, robust = robust)
  list(result = extract_extended(fit0, fit, tr, j), fit0 = fit0, fit = fit, treat = tr,
       design = design, coef_index = j)
}

fit_group_aware <- function(mat, md, treat_fc = fc_main) {
  design <- build_design(md, include_run = TRUE)
  j <- match("tissue_groupcancer", colnames(design))
  group <- factor(md$tissue_group, levels = c("normal", "cancer"))
  weighted <- voomaByGroup(mat, group = group, design = design, correlation = 0)
  fit0 <- lmFit(weighted, design)
  fit <- eBayes(fit0, trend = FALSE, robust = TRUE)
  tr <- treat(fit0, lfc = log2(treat_fc), trend = FALSE, robust = TRUE)
  list(result = extract_extended(fit0, fit, tr, j), fit0 = fit0, fit = fit, treat = tr,
       design = design, coef_index = j, weight_range = range(weighted$weights))
}

annotate_gene_result <- function(dt) {
  setnames(copy(dt), "feature_id", "entrez_id")[gene_ann, on = "entrez_id"]
}

# Mean-variance diagnostic is outcome-aware only with respect to the predeclared groups;
# it does not select genes or samples. Group-aware weights are run as a technical sensitivity.
group_stats <- rbindlist(lapply(c("normal", "cancer"), function(g) {
  idx <- metadata$tissue_group == g
  data.table(entrez_id = rownames(expr_rep), tissue_group = g,
             group_mean = rowMeans(expr_rep[, idx, drop = FALSE]),
             group_variance = apply(expr_rep[, idx, drop = FALSE], 1, var))
}))
group_stats[, log2_variance := log2(pmax(group_variance, 1e-8))]
group_stats[, mean_bin := pmin(30L, ceiling(30 * frank(group_mean, ties.method = "average") / .N)),
            by = tissue_group]
mv_binned <- group_stats[!is.na(mean_bin), .(mean_expression = median(group_mean),
                                              median_log2_variance = median(log2_variance),
                                              genes_in_bin = .N),
                          by = .(tissue_group, mean_bin)]
mv_summary <- group_stats[, .(
  genes = .N,
  spearman_mean_variance = cor(group_mean, log2_variance, method = "spearman"),
  linear_slope_log2_variance_per_expression_unit = coef(lm(log2_variance ~ group_mean))[2],
  median_log2_variance = median(log2_variance)
), by = tissue_group]
write_tsv(group_stats, "results/tables/mean_variance_gene_diagnostic.tsv")
write_tsv(mv_summary, "results/tables/mean_variance_summary.tsv")
write_tsv(mv_binned, "results/tables/figure_S01_source.tsv")
f_s01 <- ggplot(mv_binned, aes(mean_expression, median_log2_variance, color = tissue_group)) +
  geom_line(aes(linetype = tissue_group), linewidth = 1.05, lineend = "round") +
  geom_point(aes(shape = tissue_group), size = 1.65) +
  scale_color_manual(values = group_colors) + scale_linetype_manual(values = c(normal = "22", cancer = "solid")) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "Group-specific mean-variance diagnostic", subtitle = "Thirty equal-frequency expression bins; median gene variance within bin",
       x = "Mean log2 expression within group", y = "Median log2 variance",
       color = "Tissue", linetype = "Tissue", shape = "Tissue") +
  theme_portfolio() + theme(legend.position = "top")
save_figure(f_s01, "S01_mean_variance_diagnostic", 9, 6)
write_caption("S01_mean_variance_diagnostic", "Observed group-specific variance trends support trend-aware primary moderation and motivate a group-aware vooma sensitivity model.")

# Locked models ----------------------------------------------------------------
message("[03] Fitting P1, S1, S2, no-trend, group-aware and median-collapse models")
p1 <- fit_standard(expr_rep, metadata, include_run = TRUE, trend = TRUE)
s1 <- fit_standard(expr_rep, metadata, include_run = FALSE, trend = TRUE)

run_counts <- metadata[, .(normal_n = sum(tissue_group == "normal"), cancer_n = sum(tissue_group == "cancer")),
                       by = processing_run_proxy]
mixed_runs <- run_counts[normal_n > 0 & cancer_n > 0, processing_run_proxy]
s2_idx <- metadata$processing_run_proxy %in% mixed_runs
stopifnot(sum(s2_idx) == 62L, length(mixed_runs) == 6L)
s2 <- fit_standard(expr_rep[, s2_idx, drop = FALSE], metadata[s2_idx], include_run = TRUE, trend = TRUE)
s4_no_trend <- fit_standard(expr_rep, metadata, include_run = TRUE, trend = FALSE)
s4_group_aware <- fit_group_aware(expr_rep, metadata)
median_model <- fit_standard(expr_med, metadata, include_run = TRUE, trend = TRUE)

p1_dt <- annotate_gene_result(p1$result)
s1_dt <- setnames(copy(s1$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                    c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                    c("entrez_id", "s1_log2FC", "s1_standard_error", "s1_bh_fdr", "s1_treat_bh_fdr"))
s2_dt <- setnames(copy(s2$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                    c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                    c("entrez_id", "s2_log2FC", "s2_standard_error", "s2_bh_fdr", "s2_treat_bh_fdr"))
s4_dt <- setnames(copy(s4_no_trend$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                    c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                    c("entrez_id", "s4_no_trend_log2FC", "s4_no_trend_standard_error", "s4_no_trend_bh_fdr", "s4_no_trend_treat_bh_fdr"))
s4g_dt <- setnames(copy(s4_group_aware$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                     c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                     c("entrez_id", "s4_group_aware_log2FC", "s4_group_aware_standard_error", "s4_group_aware_bh_fdr", "s4_group_aware_treat_bh_fdr"))
med_dt <- setnames(copy(median_model$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                   c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                   c("entrez_id", "median_collapse_log2FC", "median_collapse_standard_error", "median_collapse_bh_fdr", "median_collapse_treat_bh_fdr"))

# No sample met the spec's signed technical-failure definition, so S3 is N/A.
# A separate, explicitly labelled influence model removes the four QC candidates.
confirmed_failures <- ledger[primary_analysis_included == FALSE, gsm_id]
qc_candidates <- ledger[decision == "retain_primary_exclude_in_qc_candidate_sensitivity", gsm_id]
stopifnot(length(confirmed_failures) == 0L, length(qc_candidates) == 4L)
sqc_idx <- !metadata$gsm_id %in% qc_candidates
sqc <- fit_standard(expr_rep[, sqc_idx, drop = FALSE], metadata[sqc_idx], include_run = TRUE, trend = TRUE)
sqc_dt <- setnames(copy(sqc$result[, .(feature_id, log2FC, standard_error, bh_fdr, treat_bh_fdr)]),
                     c("feature_id", "log2FC", "standard_error", "bh_fdr", "treat_bh_fdr"),
                     c("entrez_id", "qc_candidate_exclusion_log2FC", "qc_candidate_exclusion_standard_error",
                       "qc_candidate_exclusion_bh_fdr", "qc_candidate_exclusion_treat_bh_fdr"))

de <- Reduce(function(x, y) merge(x, y, by = "entrez_id", all.x = TRUE, sort = FALSE),
             list(p1_dt, s1_dt, s2_dt, s4_dt, s4g_dt, med_dt, sqc_dt))
de[, `:=`(s3_status = "not_applicable_no_confirmed_technical_failures",
          s3_log2FC = NA_real_, s3_standard_error = NA_real_, s3_bh_fdr = NA_real_, s3_treat_bh_fdr = NA_real_)]

# Pre-specified 1.1- and 1.5-fold TREAT sensitivities use the exact P1 fit/design.
treat_counts <- rbindlist(lapply(sort(unique(c(fc_sens, fc_main))), function(fc) {
  tr <- treat(p1$fit0, lfc = log2(fc), trend = TRUE, robust = TRUE)
  fdr <- p.adjust(tr$p.value[, p1$coef_index], method = "BH")
  data.table(minimum_fold_change = fc, minimum_abs_log2FC = log2(fc),
             tested_genes = length(fdr), treat_fdr_lt_0_05 = sum(fdr < bh_cut),
             up = sum(fdr < bh_cut & p1$result$log2FC > 0),
             down = sum(fdr < bh_cut & p1$result$log2FC < 0))
}))
write_tsv(treat_counts, "results/tables/treat_threshold_sensitivity_counts.tsv")

# Probe-level technical supplement and within-gene direction audit.
message("[03] Fitting eligible probes and auditing constituent-probe direction")
probe_ids <- pairs$probe_id
expr_probe <- expr_probe_all[probe_ids, , drop = FALSE]
probe_model <- fit_standard(expr_probe, metadata, include_run = TRUE, trend = TRUE)
probe_dt <- setnames(copy(probe_model$result), "feature_id", "probe_id")
probe_dt <- pairs[probe_dt, on = "probe_id"]
rep_sign <- de[, .(entrez_id, representative_log2FC = log2FC, representative_direction = sign(log2FC))]
probe_dt <- rep_sign[probe_dt, on = "entrez_id"]
probe_dt[, `:=`(
  same_direction_as_representative = sign(log2FC) == representative_direction,
  material_opposite_treat_signal = sign(log2FC) == -representative_direction & treat_bh_fdr < bh_cut
)]
probe_concordance <- probe_dt[, .(
  constituent_probe_count_modelled = .N,
  constituent_same_direction_count = sum(same_direction_as_representative),
  constituent_same_direction_fraction = mean(same_direction_as_representative),
  any_opposite_direction_probe = any(!same_direction_as_representative),
  any_material_opposite_treat_probe = any(material_opposite_treat_signal),
  material_opposite_probe_ids = paste(sort(probe_id[material_opposite_treat_signal]), collapse = ";")
), by = entrez_id]
fwrite(probe_dt, "results/tables/de_probe_level_supplement.tsv.gz", sep = "\t", na = "NA", quote = FALSE, compress = "gzip")

de <- merge(de, probe_concordance, by = "entrez_id", all.x = TRUE, sort = FALSE)
de[, `:=`(
  statistically_discovered = bh_fdr < bh_cut,
  priority_effect = treat_bh_fdr < bh_cut,
  representative_median_direction_concordant = sign(log2FC) == sign(median_collapse_log2FC),
  core_sensitivity_direction_concordant = sign(log2FC) == sign(s1_log2FC) &
    sign(log2FC) == sign(s2_log2FC) & sign(log2FC) == sign(s4_no_trend_log2FC) &
    sign(log2FC) == sign(s4_group_aware_log2FC) & sign(log2FC) == sign(qc_candidate_exclusion_log2FC)
)]
de[, headline_eligible := priority_effect & representative_median_direction_concordant &
     core_sensitivity_direction_concordant & !any_material_opposite_treat_probe]
de[, contrast := "cancer-minus-normal"]
de[, result_annotation_status := "current_annotation_at_run"]
de <- de[order(treat_bh_fdr, p_value, -abs(log2FC), entrez_id)]
write_tsv(de, "results/tables/de_primary_all_genes.tsv")
write_tsv(de[priority_effect == TRUE], "results/tables/de_priority_genes.tsv")
write_tsv(probe_concordance, "results/tables/probe_direction_concordance.tsv")

# Stability summaries ----------------------------------------------------------
stability_row <- function(name, other_dt, other_col, n_samples, purpose, applicable = TRUE) {
  if (!applicable) return(data.table(
    comparison = name, status = "not_applicable", purpose = purpose, n_samples = NA_integer_,
    spearman_all_tested = NA_real_, pearson_all_tested = NA_real_, sign_concordance_all_tested = NA_real_,
    p1_priority_sign_concordance = NA_real_, median_abs_effect_ratio_p1_priority = NA_real_,
    p1_priority_reversed_count = NA_integer_))
  x <- de$log2FC
  y <- other_dt[[other_col]][match(de$entrez_id, other_dt$entrez_id)]
  pri <- de$priority_effect
  data.table(
    comparison = name, status = "completed", purpose = purpose, n_samples = n_samples,
    spearman_all_tested = cor(x, y, method = "spearman"),
    pearson_all_tested = cor(x, y, method = "pearson"),
    sign_concordance_all_tested = mean(sign(x) == sign(y)),
    p1_priority_sign_concordance = mean(sign(x[pri]) == sign(y[pri])),
    median_abs_effect_ratio_p1_priority = median(abs(y[pri]) / pmax(abs(x[pri]), 1e-12)),
    p1_priority_reversed_count = sum(sign(x[pri]) != sign(y[pri]))
  )
}

stability <- rbindlist(list(
  stability_row("P1_vs_S1_group_only", s1_dt, "s1_log2FC", nrow(metadata), "run-adjustment influence"),
  stability_row("P1_vs_S2_mixed_runs", s2_dt, "s2_log2FC", sum(s2_idx), "within-run overlap evidence"),
  stability_row("P1_vs_S3_confirmed_failure_exclusion", NULL, NULL, NA_integer_, "signed technical-failure exclusion", FALSE),
  stability_row("P1_vs_QC_candidate_exclusion", sqc_dt, "qc_candidate_exclusion_log2FC", sum(sqc_idx), "candidate influence only"),
  stability_row("P1_vs_S4_no_trend", s4_dt, "s4_no_trend_log2FC", nrow(metadata), "moderation trend choice"),
  stability_row("P1_vs_S4_group_aware", s4g_dt, "s4_group_aware_log2FC", nrow(metadata), "heteroskedastic group-aware weights"),
  stability_row("P1_vs_median_collapse", med_dt, "median_collapse_log2FC", nrow(metadata), "duplicate-probe summarization")
), fill = TRUE)

green_r <- as.numeric(cfg$stability$green_spearman_min)
green_s <- as.numeric(cfg$stability$green_sign_concordance_min)
amber_r <- as.numeric(cfg$stability$amber_spearman_min)
amber_s <- as.numeric(cfg$stability$amber_sign_concordance_min)
gate_rows <- stability[comparison %in% c("P1_vs_S1_group_only", "P1_vs_S2_mixed_runs")]
traffic <- if (all(gate_rows$spearman_all_tested >= green_r & gate_rows$p1_priority_sign_concordance >= green_s)) {
  "Green"
} else if (any(gate_rows$spearman_all_tested < amber_r | gate_rows$p1_priority_sign_concordance < amber_s)) {
  "Red"
} else {
  "Amber"
}
stability[, design_traffic_light_before_leave_one_run_out := traffic]

# Leave-one-run-out influence uses parallel fork workers and returns effects only.
message("[03] Running 13 leave-one-run-out fits")
run_levels <- sort(unique(metadata$processing_run_proxy))
worker_n <- min(8L, parallel::detectCores(), length(run_levels))
loro <- mclapply(run_levels, function(run_drop) {
  keep <- metadata$processing_run_proxy != run_drop
  md <- metadata[keep]
  result <- tryCatch(fit_standard(expr_rep[, keep, drop = FALSE], md, include_run = TRUE, trend = TRUE)$result,
                     error = function(e) e)
  if (inherits(result, "error")) return(list(run = run_drop, error = conditionMessage(result)))
  list(run = run_drop, effect = setNames(result$log2FC, result$feature_id), n = sum(keep), error = NA_character_)
}, mc.cores = worker_n, mc.preschedule = FALSE)

priority_ids <- de[priority_effect == TRUE, entrez_id]
loro_rows <- rbindlist(lapply(loro, function(z) {
  rc <- run_counts[processing_run_proxy == z$run]
  if (!is.na(z$error)) return(data.table(removed_run = z$run, status = "failed", error = z$error))
  y <- z$effect[de$entrez_id]
  pri <- de$entrez_id %in% priority_ids
  data.table(
    removed_run = z$run, status = "completed", error = NA_character_, samples_retained = z$n,
    removed_normal_n = rc$normal_n, removed_cancer_n = rc$cancer_n,
    removed_run_type = ifelse(rc$normal_n > 0 & rc$cancer_n > 0, "mixed", "cancer_only"),
    spearman_all_tested = cor(de$log2FC, y, method = "spearman"),
    pearson_all_tested = cor(de$log2FC, y),
    p1_priority_sign_concordance = mean(sign(de$log2FC[pri]) == sign(y[pri])),
    p1_priority_reversed_count = sum(sign(de$log2FC[pri]) != sign(y[pri])),
    median_abs_effect_ratio_p1_priority = median(abs(y[pri]) / pmax(abs(de$log2FC[pri]), 1e-12))
  )
}), fill = TRUE)
stopifnot(nrow(loro_rows) == length(run_levels), all(loro_rows$status == "completed"))
loro_effect <- do.call(cbind, lapply(loro, function(z) z$effect[de$entrez_id]))
colnames(loro_effect) <- run_levels
de[, leave_one_run_all_same_direction := rowSums(sign(loro_effect) == sign(log2FC)) == ncol(loro_effect)]
de[, leave_one_run_min_abs_log2FC := apply(abs(loro_effect), 1, min)]

loro_downgrade <- any(loro_rows[removed_run_type == "mixed", p1_priority_sign_concordance < green_s])
traffic_final <- traffic
if (traffic_final == "Green" && loro_downgrade) traffic_final <- "Amber"
if (any(loro_rows$spearman_all_tested < amber_r | loro_rows$p1_priority_sign_concordance < amber_s)) traffic_final <- "Red"
stability[, final_design_traffic_light := traffic_final]
write_tsv(stability, "results/tables/de_sensitivity_summary.tsv")
write_tsv(loro_rows, "results/tables/leave_one_run_out_summary.tsv")

# Re-write DE tables with leave-one-run-out audit fields.
de <- de[order(treat_bh_fdr, p_value, -abs(log2FC), entrez_id)]
write_tsv(de, "results/tables/de_primary_all_genes.tsv")
write_tsv(de[priority_effect == TRUE], "results/tables/de_priority_genes.tsv")

analysis_sets <- data.table(
  analysis_id = c("P1", "S1", "S2", "S3", "SQC", "S4a", "S4b", "MDC", "PROBE", "LORO"),
  status = c(rep("completed", 3), "not_applicable", rep("completed", 6)),
  samples = c(121L, 121L, 62L, NA_integer_, 117L, 121L, 121L, 121L, 121L, NA_integer_),
  model_or_scope = c(
    "run proxy + tissue group; trend+robust", "tissue group only; trend+robust",
    "six mixed runs; run proxy + tissue group; trend+robust",
    "no signed technical failures", "exclude four reviewed QC candidates; influence only",
    "P1 with trend=FALSE; robust=TRUE", "P1 with group-specific vooma weights; robust=TRUE",
    "median collapse; P1 model", "eligible unique-Entrez probes; P1 model", "remove each run and refit P1"
  ),
  inferential_role = c("primary", rep("sensitivity", 7), "technical supplement", "influence diagnostic")
)
write_tsv(analysis_sets, "results/tables/de_analysis_set_manifest.tsv")

model_diag <- data.table(
  metric = c("tested_genes", "primary_samples", "normal_samples", "cancer_samples", "primary_design_columns",
             "primary_design_rank", "primary_residual_df", "mixed_runs", "mixed_run_samples",
             "priority_genes", "statistically_discovered_genes", "headline_eligible_genes",
             "group_aware_weight_min", "group_aware_weight_max", "final_design_traffic_light"),
  value = c(nrow(de), nrow(metadata), sum(metadata$tissue_group == "normal"), sum(metadata$tissue_group == "cancer"),
            ncol(p1$design), qr(p1$design)$rank, unique(p1$fit$df.residual)[1], length(mixed_runs), sum(s2_idx),
            sum(de$priority_effect), sum(de$statistically_discovered), sum(de$headline_eligible),
            s4_group_aware$weight_range[1], s4_group_aware$weight_range[2], traffic_final)
)
write_tsv(model_diag, "results/tables/de_model_diagnostics.tsv")

# Core figures -----------------------------------------------------------------
stability_plot <- rbindlist(list(
  de[, .(entrez_id, symbol_at_run, p1_log2FC = log2FC, sensitivity_log2FC = s1_log2FC,
         comparison = "S1: group only", priority_effect)],
  de[, .(entrez_id, symbol_at_run, p1_log2FC = log2FC, sensitivity_log2FC = s2_log2FC,
         comparison = "S2: mixed runs", priority_effect)],
  de[, .(entrez_id, symbol_at_run, p1_log2FC = log2FC, sensitivity_log2FC = qc_candidate_exclusion_log2FC,
         comparison = "QC-candidate exclusion", priority_effect)]
))
stability_labels <- stability[comparison %in% c("P1_vs_S1_group_only", "P1_vs_S2_mixed_runs", "P1_vs_QC_candidate_exclusion"),
                              .(comparison = c("S1: group only", "S2: mixed runs", "QC-candidate exclusion"),
                                label = sprintf("rho = %.3f; priority sign = %.1f%%", spearman_all_tested,
                                                100 * p1_priority_sign_concordance))]
f08 <- ggplot(stability_plot, aes(p1_log2FC, sensitivity_log2FC, color = priority_effect)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = portfolio_colors[["neutral_light"]]) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = portfolio_colors[["neutral_light"]]) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.7,
              color = portfolio_colors[["stable_dark"]]) +
  geom_point(data = stability_plot[priority_effect == FALSE], color = portfolio_colors[["neutral_mid"]],
             size = 0.45, alpha = 0.14) +
  geom_point(data = stability_plot[priority_effect == TRUE], color = portfolio_colors[["cancer"]],
             size = 0.65, alpha = 0.42) + facet_wrap(~ comparison, nrow = 1) +
  geom_text(data = stability_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.3, size = 3.0, color = portfolio_colors[["neutral_dark"]]) +
  labs(title = paste0("Design and sample-influence stability: ", traffic_final),
       x = "P1 adjusted log2 fold-change", y = "Sensitivity log2 fold-change") +
  coord_equal() + theme_portfolio() + theme(legend.position = "none")
save_figure(f08, "F08_design_stability", 14, 5.8)
write_tsv(stability_plot, "results/tables/figure_F08_source.tsv")
write_caption("F08_design_stability", "Each panel compares cancer-minus-normal effects with the locked P1 model. Sign concordance is evaluated on P1 TREAT-priority genes; S2 significance retention is not required.")

f09_data <- de[, .(entrez_id, symbol_at_run, average_expression, log2FC, bh_fdr, treat_bh_fdr,
                   result_class = fifelse(priority_effect, "TREAT priority",
                                          fifelse(statistically_discovered, "BH discovery", "Not discovered")))]
f09_data[, visual_class := fcase(
  result_class == "TREAT priority" & log2FC > 0, "Priority: cancer-up",
  result_class == "TREAT priority" & log2FC < 0, "Priority: normal-up",
  result_class == "BH discovery", "BH discovery",
  default = "Not discovered"
)]
result_colors <- c(
  "Not discovered" = portfolio_colors[["neutral_light"]],
  "BH discovery" = portfolio_colors[["stable"]],
  "Priority: normal-up" = portfolio_colors[["normal"]],
  "Priority: cancer-up" = portfolio_colors[["cancer"]]
)
f09 <- ggplot(f09_data, aes(average_expression, log2FC, color = visual_class)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = portfolio_colors[["neutral_mid"]]) +
  geom_point(data = f09_data[visual_class == "Not discovered"], size = 0.55, alpha = 0.32) +
  geom_point(data = f09_data[visual_class == "BH discovery"], size = 0.65, alpha = 0.46) +
  geom_point(data = f09_data[grepl("^Priority", visual_class)], size = 0.78, alpha = 0.62) +
  scale_color_manual(values = result_colors, breaks = names(result_colors)) +
  labs(title = "P1 mean-difference plot", x = "Average log2 expression", y = "Adjusted cancer-minus-normal log2FC", color = NULL) +
  theme_portfolio() + theme(legend.position = "top")
save_figure(f09, "F09_MA_plot", 9, 6.5)
write_tsv(f09_data, "results/tables/figure_F09_source.tsv")
write_caption("F09_MA_plot", "P1 effects are adjusted for the processing-run proxy. TREAT priority directly tests whether the absolute effect exceeds log2(1.2).")

f10_data <- de[, .(entrez_id, symbol_at_run, log2FC, p_value, bh_fdr, treat_bh_fdr,
                   result_class = fifelse(priority_effect, "TREAT priority",
                                          fifelse(statistically_discovered, "BH discovery", "Not discovered")))]
f10_data[, minus_log10_p := -log10(pmax(p_value, .Machine$double.xmin))]
f10_data[, visual_class := fcase(
  result_class == "TREAT priority" & log2FC > 0, "Priority: cancer-up",
  result_class == "TREAT priority" & log2FC < 0, "Priority: normal-up",
  result_class == "BH discovery", "BH discovery",
  default = "Not discovered"
)]
label_dt <- rbind(
  f10_data[result_class == "TREAT priority" & log2FC > 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, cfg$reporting$volcano_labels_each_direction)],
  f10_data[result_class == "TREAT priority" & log2FC < 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, cfg$reporting$volcano_labels_each_direction)]
)
f10 <- ggplot(f10_data, aes(log2FC, minus_log10_p, color = visual_class)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = portfolio_colors[["neutral_light"]]) +
  geom_vline(xintercept = c(-lfc_main, lfc_main), linetype = "22", linewidth = 0.55,
             color = portfolio_colors[["neutral_mid"]]) +
  geom_point(data = f10_data[visual_class == "Not discovered"], size = 0.55, alpha = 0.30) +
  geom_point(data = f10_data[visual_class == "BH discovery"], size = 0.65, alpha = 0.46) +
  geom_point(data = f10_data[grepl("^Priority", visual_class)], size = 0.8, alpha = 0.64) +
  ggrepel::geom_text_repel(data = label_dt, aes(label = symbol_at_run), size = 2.6,
                           fontface = "italic", color = portfolio_colors[["neutral_dark"]],
                           segment.color = portfolio_colors[["neutral_mid"]], segment.size = 0.28,
                           box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf, seed = cfg$project$seed,
                           min.segment.length = 0) +
  scale_color_manual(values = result_colors, breaks = names(result_colors)) +
  labs(title = "P1 differential-expression evidence", subtitle = "Dashed lines show the 1.2-fold minimum effect; classes use BH-adjusted tests",
       x = "Adjusted cancer-minus-normal log2FC", y = expression(-log[10](raw~P)), color = NULL) +
  theme_portfolio() + theme(legend.position = "top")
save_figure(f10, "F10_volcano_plot", 10, 7)
write_tsv(f10_data, "results/tables/figure_F10_source.tsv")
write_caption("F10_volcano_plot", "The y-axis uses raw moderated-test P values for display, while all discovery and priority labels use BH FDR. The TREAT class tests an absolute effect greater than 1.2-fold.")

heat_candidates <- de[headline_eligible == TRUE]
if (!nrow(heat_candidates)) heat_candidates <- de[priority_effect == TRUE]
up_ids <- heat_candidates[log2FC > 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, cfg$reporting$heatmap_each_direction), entrez_id]
down_ids <- heat_candidates[log2FC < 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, cfg$reporting$heatmap_each_direction), entrez_id]
heat_ids <- c(up_ids, down_ids)
heat_raw <- expr_rep[heat_ids, , drop = FALSE]
heat_z <- t(scale(t(heat_raw)))
heat_z[heat_z > 2.5] <- 2.5; heat_z[heat_z < -2.5] <- -2.5
heat_labels <- gene_ann$symbol_at_run[match(heat_ids, gene_ann$entrez_id)]
heat_labels[is.na(heat_labels) | heat_labels == ""] <- heat_ids[is.na(heat_labels) | heat_labels == ""]
rownames(heat_z) <- make.unique(heat_labels)
ann_col <- data.frame(Tissue = factor(metadata$tissue_group, levels = c("normal", "cancer")),
                      Run = factor(metadata$processing_run_proxy), row.names = metadata$gsm_id)
run_colors <- portfolio_run_colors(levels(ann_col$Run))
ph <- pheatmap(heat_z, color = portfolio_diverging(101),
               breaks = seq(-2.5, 2.5, length.out = 102), cluster_rows = TRUE, cluster_cols = TRUE,
               clustering_method = "average", annotation_col = ann_col,
               annotation_colors = list(Tissue = group_colors, Run = run_colors),
               show_colnames = FALSE, border_color = NA, main = "Priority-gene expression across samples",
               fontsize = 8.5, fontsize_row = 8, silent = TRUE)
draw_f11 <- function() { grid::grid.newpage(); grid::grid.draw(ph$gtable) }
save_base_figure("F11_priority_gene_heatmap", draw_f11, 13, 10)
heat_source <- as.data.table(heat_raw, keep.rownames = "entrez_id")
write_tsv(heat_source, "results/tables/figure_F11_source.tsv")
write_caption("F11_priority_gene_heatmap", "Rows are selected from headline-eligible TREAT-priority genes without using the heatmap for inference. Values are row-wise z-scores capped at +/-2.5; source data preserve original log2 expression.")

forest_base <- de[headline_eligible == TRUE]
if (!nrow(forest_base)) forest_base <- de[priority_effect == TRUE]
forest_ids <- c(forest_base[log2FC > 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, 10), entrez_id],
                forest_base[log2FC < 0][order(treat_bh_fdr, -abs(log2FC))][1:min(.N, 10), entrez_id])
forest <- de[entrez_id %in% forest_ids]
forest[, label := fifelse(is.na(symbol_at_run) | symbol_at_run == "", entrez_id, symbol_at_run)]
forest[, direction := fifelse(log2FC > 0, "Cancer-up", "Normal-up")]
forest[, label := factor(label, levels = label[order(log2FC)])]
forest_points <- melt(forest[, .(entrez_id, label, P1 = log2FC, S1 = s1_log2FC, S2 = s2_log2FC,
                                 `QC candidate exclusion` = qc_candidate_exclusion_log2FC,
                                 `S4 no trend` = s4_no_trend_log2FC)],
                      id.vars = c("entrez_id", "label"), variable.name = "analysis", value.name = "effect")
f12 <- ggplot(forest, aes(log2FC, label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.55, color = portfolio_colors[["neutral_mid"]]) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high), orientation = "y", width = 0.20,
                linewidth = 0.55, color = portfolio_colors[["neutral_dark"]]) +
  geom_point(aes(fill = direction), shape = 21, size = 2.8, stroke = 0.55,
             color = portfolio_colors[["neutral_dark"]]) +
  geom_point(data = forest_points[analysis != "P1"], aes(effect, label, color = analysis, shape = analysis),
             inherit.aes = FALSE, size = 1.85, alpha = 0.9) +
  scale_fill_manual(values = c("Cancer-up" = portfolio_colors[["cancer"]],
                               "Normal-up" = portfolio_colors[["normal"]]), name = "P1 direction") +
  scale_color_manual(values = c(
    S1 = portfolio_colors[["cancer_dark"]],
    S2 = portfolio_colors[["normal_dark"]],
    `QC candidate exclusion` = portfolio_colors[["stable_dark"]],
    `S4 no trend` = portfolio_colors[["neutral_mid"]]
  )) +
  scale_shape_manual(values = c(S1 = 16, S2 = 17, `QC candidate exclusion` = 15, `S4 no trend` = 18)) +
  labs(title = "Headline effect estimates and sensitivity", x = "Cancer-minus-normal log2FC", y = NULL,
       color = "Sensitivity", shape = "Sensitivity") +
  theme_portfolio() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_figure(f12, "F12_top_effect_forest", 10, 8)
write_tsv(merge(forest[, .(entrez_id, label, log2FC, standard_error, ci95_low, ci95_high)], forest_points,
                by = c("entrez_id", "label"), all.x = TRUE), "results/tables/figure_F12_source.tsv")
write_caption("F12_top_effect_forest", "Black intervals are P1 95% confidence intervals; colored points show pre-specified sensitivity effects. Genes are selected from headline-eligible TREAT-priority results.")

message(sprintf("[03] PASS: %d tested; %d BH discoveries; %d TREAT priorities; %d headline eligible; stability %s",
                nrow(de), sum(de$statistically_discovered), sum(de$priority_effect), sum(de$headline_eligible), traffic_final))

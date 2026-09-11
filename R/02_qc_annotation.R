source("R/common.R")
ctx <- init_project()
cfg <- ctx$cfg
ensure_dirs()

suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(MASS)
  library(parallel)
})

message("[02] Running processed-data QC and locked probe annotation")
expr <- readRDS("data/derived/expression_probe_matrix.rds")
metadata <- as.data.table(readRDS("data/derived/sample_metadata.rds"))
stopifnot(identical(colnames(expr), metadata$gsm_id))

# Probe-level annotation audit -------------------------------------------------
probe_ids <- rownames(expr)
probe_variance <- apply(expr, 1, var)
probe_median <- apply(expr, 1, median)
constant <- probe_variance == 0

ann <- suppressMessages(AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("ENTREZID"),
  keytype = "PROBEID"
))
ann <- unique(as.data.table(ann)[, .(probe_id = PROBEID, entrez_id = as.character(ENTREZID))])
map_pairs <- unique(ann[!is.na(entrez_id) & nzchar(entrez_id)])
map_count <- map_pairs[, .(entrez_count = uniqueN(entrez_id),
                           mapped_entrez_ids = paste(sort(unique(entrez_id)), collapse = ";")), by = probe_id]

probe_audit <- data.table(
  probe_id = probe_ids,
  probe_variance = probe_variance,
  pooled_median_expression = probe_median,
  is_constant = constant,
  is_affx_control = startsWith(probe_ids, "AFFX-")
)
probe_audit[, suffix_class := fcase(
  is_affx_control, "control",
  grepl("_x_at$", probe_id), "x_at",
  grepl("_s_at$", probe_id), "s_at",
  grepl("_a_at$", probe_id), "a_at",
  grepl("_at$", probe_id), "at",
  default = "unknown"
)]
probe_audit <- map_count[probe_audit, on = "probe_id"]
probe_audit[is.na(entrez_count), entrez_count := 0L]
probe_audit[is.na(mapped_entrez_ids), mapped_entrez_ids := NA_character_]
probe_audit[, mapping_status := fcase(
  is_affx_control, "control_affx",
  is_constant, "constant",
  entrez_count == 0L, "unmapped",
  entrez_count > 1L, "one_to_many",
  suffix_class == "x_at", "cross_hybridizing_x_at",
  suffix_class == "unknown", "unknown_suffix",
  default = "primary_eligible"
)]
probe_audit[, primary_eligible := mapping_status == "primary_eligible"]

eligible_denominator <- probe_audit[!is_affx_control & !is_constant, .N]
eligible_fraction <- probe_audit[primary_eligible == TRUE, .N] / eligible_denominator
if (eligible_fraction < cfg$annotation$minimum_eligible_probe_fraction) {
  stop(sprintf("Hard stop: eligible probe fraction %.3f is below %.3f", eligible_fraction,
               cfg$annotation$minimum_eligible_probe_fraction))
}

eligible_pairs <- map_pairs[probe_audit[primary_eligible == TRUE, .(probe_id)], on = "probe_id", nomatch = 0]
eligible_pairs <- unique(eligible_pairs)
eligible_pairs <- probe_audit[, .(probe_id, suffix_class, pooled_median_expression)][eligible_pairs, on = "probe_id"]
eligible_pairs[, specificity_tier := fcase(suffix_class == "at", 1L,
                                            suffix_class == "a_at", 2L,
                                            suffix_class == "s_at", 3L,
                                            default = 99L)]
setorder(eligible_pairs, entrez_id, specificity_tier, -pooled_median_expression, probe_id)
eligible_pairs[, selection_rank := seq_len(.N), by = entrez_id]
eligible_pairs[, selected_representative := selection_rank == 1L]
eligible_pairs[, eligible_probe_count_for_gene := .N, by = entrez_id]

representatives <- eligible_pairs[selected_representative == TRUE]
if (uniqueN(representatives$entrez_id) < cfg$annotation$minimum_unique_entrez) {
  stop(sprintf("Hard stop: only %d unique Entrez genes are testable", uniqueN(representatives$entrez_id)))
}

gene_ids <- representatives$entrez_id
gene_info <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENTREZID"
))
gene_info <- as.data.table(gene_info)
gene_info[, ENTREZID := as.character(ENTREZID)]
setorder(gene_info, ENTREZID, SYMBOL, GENENAME, na.last = TRUE)
gene_info <- gene_info[, .(
  symbol_at_run = {x <- unique(na.omit(SYMBOL)); if (length(x)) x[[1]] else NA_character_},
  gene_name_at_run = {x <- unique(na.omit(GENENAME)); if (length(x)) x[[1]] else NA_character_}
), by = .(entrez_id = ENTREZID)]

rep_indices <- match(representatives$probe_id, rownames(expr))
gene_matrix <- expr[rep_indices, , drop = FALSE]
rownames(gene_matrix) <- representatives$entrez_id

split_probe_ids <- split(eligible_pairs$probe_id, eligible_pairs$entrez_id)
all_gene_ids <- names(split_probe_ids)
ncores <- min(8L, parallel::detectCores(logical = FALSE), length(all_gene_ids))
chunks <- split(all_gene_ids, cut(seq_along(all_gene_ids), breaks = ncores, labels = FALSE))
median_parts <- parallel::mclapply(chunks, function(ids) {
  ans <- matrix(NA_real_, nrow = length(ids), ncol = ncol(expr),
                dimnames = list(ids, colnames(expr)))
  for (i in seq_along(ids)) {
    idx <- match(split_probe_ids[[ids[[i]]]], rownames(expr))
    ans[i, ] <- if (length(idx) == 1L) expr[idx, ] else apply(expr[idx, , drop = FALSE], 2, median)
  }
  ans
}, mc.cores = ncores)
median_gene_matrix <- do.call(rbind, median_parts)
median_gene_matrix <- median_gene_matrix[rownames(gene_matrix), , drop = FALSE]

gene_annotation <- representatives[, .(
  entrez_id, representative_probe_id = probe_id, specificity_tier, suffix_class,
  representative_pooled_median = pooled_median_expression, eligible_probe_count_for_gene
)]
gene_annotation <- gene_info[gene_annotation, on = "entrez_id"]
gene_annotation[, annotation_package := "hgu133plus2.db"]
gene_annotation[, annotation_package_version := as.character(packageVersion("hgu133plus2.db"))]
gene_annotation[, human_gene_package_version := as.character(packageVersion("org.Hs.eg.db"))]

probe_audit[, many_probes_per_gene := FALSE]
multi_gene_ids <- eligible_pairs[, .N, by = entrez_id][N > 1L, entrez_id]
probe_audit[probe_id %in% eligible_pairs[entrez_id %in% multi_gene_ids, probe_id], many_probes_per_gene := TRUE]

mapping_summary <- rbind(
  probe_audit[, .(count = .N), by = .(metric = paste0("mapping_status:", mapping_status))],
  probe_audit[, .(count = .N), by = .(metric = paste0("suffix_class:", suffix_class))],
  data.table(metric = c("total_probe_sets", "noncontrol_nonconstant_denominator", "primary_eligible_probes",
                        "primary_eligible_fraction", "unique_testable_entrez", "genes_with_multiple_eligible_probes"),
             count = c(nrow(probe_audit), eligible_denominator, probe_audit[primary_eligible == TRUE, .N],
                       eligible_fraction, nrow(gene_annotation), length(multi_gene_ids)))
, fill = TRUE)

write_tsv(probe_audit, "results/tables/probe_mapping_audit.tsv")
write_tsv(eligible_pairs, "results/tables/probe_selection_audit.tsv")
write_tsv(gene_annotation, "results/tables/gene_annotation.tsv")
write_tsv(mapping_summary, "results/tables/annotation_metrics.tsv")
saveRDS(gene_matrix, "data/derived/expression_gene_representative.rds", compress = "xz")
saveRDS(median_gene_matrix, "data/derived/expression_gene_median_collapse.rds", compress = "xz")
saveRDS(gene_annotation, "data/derived/gene_annotation.rds", compress = "xz")
saveRDS(eligible_pairs, "data/derived/eligible_probe_pairs.rds", compress = "xz")

# Sample-level distribution QC ------------------------------------------------
sample_metrics <- data.table(
  gsm_id = colnames(expr),
  sample_median = apply(expr, 2, median),
  sample_iqr = apply(expr, 2, IQR),
  sample_mad = apply(expr, 2, mad),
  q01 = apply(expr, 2, quantile, probs = 0.01),
  q99 = apply(expr, 2, quantile, probs = 0.99),
  dynamic_range = apply(expr, 2, function(x) diff(range(x)))
)
for (nm in c("sample_median", "sample_iqr", "sample_mad", "q01", "q99", "dynamic_range")) {
  sample_metrics[[paste0(nm, "_robust_z")]] <- robust_z(sample_metrics[[nm]])
}

gene_variability <- apply(gene_matrix, 1, mad)
top5000_ids <- names(sort(gene_variability, decreasing = TRUE))[seq_len(min(cfg$qc$pca_top_variable_genes, nrow(gene_matrix)))]
top1000_ids <- names(sort(gene_variability, decreasing = TRUE))[seq_len(min(cfg$qc$pca_sensitivity_genes, nrow(gene_matrix)))]

run_pca <- function(ids, label) {
  pca <- prcomp(t(gene_matrix[ids, , drop = FALSE]), center = TRUE, scale. = FALSE)
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  scores <- as.data.table(pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE], keep.rownames = "gsm_id")
  scores[, pca_feature_set := label]
  list(pca = pca, variance = variance, scores = scores)
}
pca5000 <- run_pca(top5000_ids, "top5000_MAD")
pca1000 <- run_pca(top1000_ids, "top1000_MAD")
pca_all <- run_pca(rownames(gene_matrix), "all_eligible_genes")

pca_scores <- rbindlist(list(pca5000$scores, pca1000$scores, pca_all$scores), fill = TRUE)
pca_variance <- rbindlist(list(
  data.table(pca_feature_set = "top5000_MAD", PC = seq_along(pca5000$variance), variance_explained = pca5000$variance),
  data.table(pca_feature_set = "top1000_MAD", PC = seq_along(pca1000$variance), variance_explained = pca1000$variance),
  data.table(pca_feature_set = "all_eligible_genes", PC = seq_along(pca_all$variance), variance_explained = pca_all$variance)
))

cumvar <- cumsum(pca5000$variance)
n_pc_distance <- min(as.integer(cfg$qc$pca_max_components_for_distance), which(cumvar >= cfg$qc$pca_cumulative_variance_for_distance)[1])
n_pc_distance <- max(2L, n_pc_distance)
pc_for_distance <- pca5000$pca$x[, seq_len(n_pc_distance), drop = FALSE]
rob <- tryCatch(MASS::cov.rob(pc_for_distance, method = "mve"), error = function(e) NULL)
if (is.null(rob)) {
  pca_distance <- mahalanobis(pc_for_distance, colMeans(pc_for_distance), cov(pc_for_distance))
  pca_distance_method <- "classical_covariance_fallback"
} else {
  pca_distance <- mahalanobis(pc_for_distance, rob$center, rob$cov)
  pca_distance_method <- "MASS_cov.rob_mve"
}
pca_cutoff <- qchisq(cfg$qc$robust_pca_probability, df = n_pc_distance)

sample_correlation <- cor(gene_matrix, method = "pearson")
sample_spearman <- cor(gene_matrix, method = "spearman")
median_corr <- apply(sample_correlation, 2, function(x) median(x[names(x) != names(which.max(x))]))
connectivity <- colSums(sample_correlation) - 1

sample_metrics[, median_pearson_correlation := median_corr[gsm_id]]
sample_metrics[, correlation_connectivity := connectivity[gsm_id]]
sample_metrics[, median_correlation_robust_z := robust_z(median_pearson_correlation)]
sample_metrics[, connectivity_robust_z := robust_z(correlation_connectivity)]
sample_metrics[, robust_pca_distance := pca_distance[gsm_id]]
sample_metrics[, robust_pca_cutoff := pca_cutoff]
sample_metrics[, pca_distance_method := pca_distance_method]

zlim <- cfg$qc$robust_z_warning
sample_metrics[, distribution_flag := apply(.SD, 1, function(x) any(abs(x) > zlim)),
               .SDcols = patterns("^(sample_median|sample_iqr|sample_mad|q01|q99|dynamic_range)_robust_z$")]
sample_metrics[, correlation_flag := median_correlation_robust_z < -zlim | connectivity_robust_z < -zlim]
sample_metrics[, pca_flag := robust_pca_distance > robust_pca_cutoff]
sample_metrics[, raw_evidence_flag := FALSE]
sample_metrics[, independent_signal_count := as.integer(distribution_flag) + as.integer(correlation_flag) +
                 as.integer(pca_flag) + as.integer(raw_evidence_flag)]
sample_metrics[, candidate_exclusion_review := independent_signal_count >= 2L]
sample_metrics <- metadata[, .(gsm_id, tissue_group, processing_run_proxy)][sample_metrics, on = "gsm_id"]

exclusion_ledger <- sample_metrics[candidate_exclusion_review == TRUE, .(
  gsm_id, tissue_group, processing_run_proxy, independent_signal_count,
  distribution_flag, correlation_flag, pca_flag, raw_evidence_flag,
  decision = "retain_pending_scientific_review",
  rationale = "No sample is automatically excluded; requires independent technical corroboration and signed review."
)]
if (!nrow(exclusion_ledger)) {
  exclusion_ledger <- data.table(
    gsm_id = NA_character_, tissue_group = NA_character_, processing_run_proxy = NA_character_,
    independent_signal_count = 0L, distribution_flag = FALSE, correlation_flag = FALSE,
    pca_flag = FALSE, raw_evidence_flag = FALSE, decision = "no_candidate_samples",
    rationale = "No sample triggered two independent processed-data QC signals."
  )
}

# Metadata-PC associations
score_primary <- as.data.table(pca5000$pca$x[, seq_len(min(10, ncol(pca5000$pca$x))), drop = FALSE], keep.rownames = "gsm_id")
score_primary <- metadata[, .(gsm_id, tissue_group, processing_run_proxy)][score_primary, on = "gsm_id"]
pc_assoc <- rbindlist(lapply(grep("^PC", names(score_primary), value = TRUE), function(pc) {
  y <- score_primary[[pc]]
  fit_group <- lm(y ~ factor(tissue_group), data = score_primary)
  fit_run <- lm(y ~ factor(processing_run_proxy), data = score_primary)
  fit_full <- lm(y ~ factor(processing_run_proxy) + factor(tissue_group), data = score_primary)
  rss0 <- sum((y - mean(y))^2)
  rss_group <- deviance(fit_group); rss_run <- deviance(fit_run); rss_full <- deviance(fit_full)
  data.table(
    PC = pc,
    group_r2 = 1 - rss_group / rss0,
    run_r2 = 1 - rss_run / rss0,
    group_partial_r2_given_run = (rss_run - rss_full) / rss_run,
    run_partial_r2_given_group = (rss_group - rss_full) / rss_group
  )
}))

# Visualization-only batch adjustment preserving tissue group
vis_design <- model.matrix(~ tissue_group, data = metadata)
adjusted_gene_matrix <- limma::removeBatchEffect(
  gene_matrix,
  batch = factor(metadata$processing_run_proxy),
  design = vis_design
)
pca_adjusted <- prcomp(t(adjusted_gene_matrix[top5000_ids, , drop = FALSE]), center = TRUE, scale. = FALSE)

write_tsv(sample_metrics, "results/tables/sample_qc_metrics.tsv")
write_tsv(exclusion_ledger, "results/tables/sample_exclusion_ledger.tsv")
write_tsv(pca_scores, "results/tables/pca_scores.tsv")
write_tsv(pca_variance, "results/tables/pca_variance.tsv")
write_tsv(pc_assoc, "results/tables/pc_metadata_associations.tsv")
write_tsv(as.data.table(sample_correlation, keep.rownames = "gsm_id"), "results/tables/sample_correlation_matrix.tsv")
write_tsv(as.data.table(sample_spearman, keep.rownames = "gsm_id"), "results/tables/sample_spearman_correlation_matrix.tsv")
saveRDS(list(top5000_ids = top5000_ids, top1000_ids = top1000_ids,
             pca5000 = pca5000$pca, pca1000 = pca1000$pca, pca_all = pca_all$pca,
             pca_adjusted = pca_adjusted, correlation = sample_correlation),
        "data/derived/qc_objects.rds", compress = "xz")

# Figures F01-F07 -------------------------------------------------------------
flow_data <- data.table(
  step = factor(c("Deposited probe sets", "Non-control, non-constant", "Eligible unambiguous probes", "Unique tested genes"),
                levels = c("Deposited probe sets", "Non-control, non-constant", "Eligible unambiguous probes", "Unique tested genes")),
  count = c(nrow(expr), eligible_denominator, probe_audit[primary_eligible == TRUE, .N], nrow(gene_matrix))
)
f01 <- ggplot(flow_data, aes(step, count, group = 1)) +
  geom_line(color = portfolio_colors[["stable_dark"]], linewidth = 1.15, lineend = "round") +
  geom_point(shape = 21, fill = portfolio_colors[["stable"]], color = "white", stroke = 1, size = 4.6) +
  geom_text(aes(label = scales::comma(count)), vjust = -1.0, fontface = "bold", size = 3.6,
            color = portfolio_colors[["neutral_dark"]]) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_discrete(labels = c(
    "Deposited probe sets" = "Deposited\nprobe sets",
    "Non-control, non-constant" = "Non-control,\nnon-constant",
    "Eligible unambiguous probes" = "Eligible unambiguous\nprobes",
    "Unique tested genes" = "Unique tested\ngenes"
  )) +
  labs(title = "Analysis feature flow", x = NULL, y = "Features") + theme_portfolio() +
  theme(axis.text.x = element_text(hjust = 0.5), axis.line.x = element_blank(), axis.ticks.x = element_blank())
save_figure(f01, "F01_analysis_feature_flow", 9, 6)
write_tsv(flow_data, "results/tables/figure_F01_source.tsv")
write_caption("F01_analysis_feature_flow", "Feature counts from deposited GPL570 probe sets to the outcome-blind representative-gene analysis matrix.")

metadata_fields <- intersect(c("age", "size", "grade", "er_status", "lymph_node_status"), names(metadata))
complete_data <- rbindlist(lapply(metadata_fields, function(field) {
  metadata[, .(complete_fraction = mean(!is.na(get(field))), n_complete = sum(!is.na(get(field))), n_total = .N),
           by = tissue_group][, field := field]
}))
complete_data[, field := factor(field, levels = rev(metadata_fields))]
f02a <- ggplot(complete_data, aes(tissue_group, field, fill = complete_fraction)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%d/%d", n_complete, n_total)), size = 3.1,
            color = portfolio_colors[["neutral_dark"]]) +
  scale_fill_gradientn(colors = c(portfolio_colors[["neutral_pale"]], portfolio_colors[["stable_light"]],
                                  portfolio_colors[["stable"]]), limits = c(0, 1), name = "Complete") +
  labs(title = "Metadata completeness", x = NULL, y = NULL) + theme_portfolio()
run_long <- melt(metadata[, .N, by = .(processing_run_proxy, tissue_group)],
                 id.vars = c("processing_run_proxy", "tissue_group"), measure.vars = "N")
f02b <- ggplot(run_long, aes(processing_run_proxy, value, fill = tissue_group)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.25) + scale_fill_manual(values = group_colors) +
  labs(title = "Tissue group by processing-run proxy", x = "Processing-run proxy", y = "Samples", fill = "Tissue") +
  theme_portfolio() + theme(axis.text.x = element_text(angle = 42, hjust = 1), legend.position = "top")
f02 <- arrangeGrob(f02a, f02b, ncol = 2, widths = c(0.8, 1.4))
save_figure(f02, "F02_metadata_and_run_structure", 13, 6)
write_tsv(complete_data, "results/tables/figure_F02_metadata_source.tsv")
write_tsv(run_long, "results/tables/figure_F02_run_source.tsv")
write_caption("F02_metadata_and_run_structure", "Normal age and tumor-only pathology fields are unavailable, while group and processing-run proxy are partially confounded but overlap in six runs.")

density_list <- lapply(seq_len(ncol(expr)), function(j) {
  d <- density(expr[, j], n = 512)
  data.table(gsm_id = colnames(expr)[j], x = d$x, density = d$y,
             tissue_group = metadata$tissue_group[j], processing_run_proxy = metadata$processing_run_proxy[j])
})
density_data <- rbindlist(density_list)
draw_f03 <- function() {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(2, 1), mar = c(4.2, 4.3, 2.5, 1.2), family = "sans",
      fg = portfolio_colors[["neutral_dark"]], col.axis = portfolio_colors[["neutral_mid"]],
      col.lab = portfolio_colors[["neutral_dark"]], las = 1, bty = "l", tcl = -0.25)
  dens_split <- split(density_data, density_data$gsm_id)
  xlim <- range(density_data$x); ylim <- range(density_data$density)
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Deposited log2 GC-RMA expression", ylab = "Density",
       main = "a  Per-sample expression densities", cex.main = 1.15, font.main = 2)
  for (id in names(dens_split)) {
    dd <- dens_split[[id]]
    lines(dd$x, dd$density, col = adjustcolor(group_colors[dd$tissue_group[1]], 0.20),
          lwd = 0.65, lty = c(normal = 2, cancer = 1)[dd$tissue_group[1]])
  }
  legend("topright", legend = c("cancer", "normal"), col = group_colors[c("cancer", "normal")],
         lwd = 2.2, lty = c(1, 2), bty = "n", cex = 0.82)
  boxplot(as.data.frame(expr), outline = FALSE, names = FALSE, xaxt = "n",
          col = unname(group_colors[metadata$tissue_group]), border = portfolio_colors[["neutral_dark"]],
          boxwex = 0.72, whisklty = 1, staplewex = 0.45, medlwd = 1.1,
          ylab = "Log2 GC-RMA expression", xlab = "Samples (matrix order)",
          main = "b  Per-sample expression distributions", cex.main = 1.15, font.main = 2)
}
save_base_figure("F03_expression_distributions", draw_f03, 14, 9)
write_tsv(density_data, "results/tables/figure_F03_density_source.tsv")
write_tsv(sample_metrics[, .(gsm_id, tissue_group, processing_run_proxy, sample_median, sample_iqr, q01, q99)],
          "results/tables/figure_F03_boxplot_source.tsv")
write_caption("F03_expression_distributions", "Deposited log2 GC-RMA values are displayed without re-normalization. Group-wide shifts are not treated as automatic technical failure.")

metric_z_cols <- grep("_robust_z$", names(sample_metrics), value = TRUE)
metric_long <- melt(sample_metrics, id.vars = c("gsm_id", "tissue_group", "processing_run_proxy"),
                    measure.vars = metric_z_cols, variable.name = "metric", value.name = "robust_z")
metric_labels <- c(
  sample_median_robust_z = "Sample median",
  sample_iqr_robust_z = "Interquartile range",
  sample_mad_robust_z = "Median absolute deviation",
  q01_robust_z = "1st percentile",
  q99_robust_z = "99th percentile",
  dynamic_range_robust_z = "Dynamic range",
  median_correlation_robust_z = "Median correlation",
  connectivity_robust_z = "Correlation connectivity"
)
f04 <- ggplot(metric_long, aes(gsm_id, robust_z, color = tissue_group)) +
  geom_hline(yintercept = c(-zlim, zlim), linetype = "22", linewidth = 0.45,
             color = portfolio_colors[["neutral_mid"]]) +
  geom_point(aes(shape = tissue_group), size = 1.25, alpha = 0.82, stroke = 0.2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2, labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = group_colors) + scale_shape_manual(values = group_shapes) +
  labs(title = "Sample-level QC metrics", x = "Samples (matrix order)", y = "Robust z-score",
       color = "Tissue", shape = "Tissue") +
  theme_portfolio() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
save_figure(f04, "F04_sample_qc_metrics", 12, 10)
write_tsv(metric_long, "results/tables/figure_F04_source.tsv")
write_caption("F04_sample_qc_metrics", "Dashed lines mark the pre-specified robust-z warning boundary. No single metric automatically excludes a sample.")

primary_scores <- as.data.table(pca5000$pca$x[, 1:3, drop = FALSE], keep.rownames = "gsm_id")
primary_scores <- metadata[, .(gsm_id, tissue_group, processing_run_proxy)][primary_scores, on = "gsm_id"]
var_pct <- 100 * pca5000$variance
run_levels <- sort(unique(primary_scores$processing_run_proxy))
run_pca_colors <- portfolio_run_colors(paste0("Run: ", run_levels))
pca_panel_data <- rbind(
  primary_scores[, .(gsm_id, PC1, PC2, panel = "a  Tissue group", display_key = paste0("Tissue: ", tissue_group))],
  primary_scores[, .(gsm_id, PC1, PC2, panel = "b  Processing-run proxy", display_key = paste0("Run: ", processing_run_proxy))]
)
pca_colors <- c(setNames(group_colors, paste0("Tissue: ", names(group_colors))), run_pca_colors)
f05 <- ggplot(pca_panel_data, aes(PC1, PC2, color = display_key)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = portfolio_colors[["neutral_light"]]) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = portfolio_colors[["neutral_light"]]) +
  geom_point(size = 2.45, alpha = 0.88) + facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(values = pca_colors) +
  labs(title = "Primary PCA: top 5,000 pooled-MAD genes", subtitle = "Genes centered; not variance-scaled",
       x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]), color = NULL) +
  theme_portfolio() + guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.2, alpha = 1)))
save_figure(f05, "F05_primary_PCA", 14, 7)
write_tsv(primary_scores, "results/tables/figure_F05_source.tsv")
write_caption("F05_primary_PCA", "Unsupervised PCA uses pooled MAD without tissue labels. The same 121 samples are shown in both panels, colored separately by tissue group and processing-run proxy.")

adjusted_scores <- as.data.table(pca_adjusted$x[, 1:3, drop = FALSE], keep.rownames = "gsm_id")
adjusted_scores <- metadata[, .(gsm_id, tissue_group, processing_run_proxy)][adjusted_scores, on = "gsm_id"]
adj_var <- 100 * pca_adjusted$sdev^2 / sum(pca_adjusted$sdev^2)
f06 <- ggplot(adjusted_scores, aes(PC1, PC2, color = tissue_group)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = portfolio_colors[["neutral_light"]]) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = portfolio_colors[["neutral_light"]]) +
  geom_point(aes(shape = tissue_group), size = 2.7, alpha = 0.88, stroke = 0.25) +
  scale_color_manual(values = group_colors) + scale_shape_manual(values = group_shapes) +
  labs(title = "Visualization-only run-adjusted PCA", subtitle = "Run removed while preserving tissue-group design; never used as the DE response",
       x = sprintf("PC1 (%.1f%%)", adj_var[1]), y = sprintf("PC2 (%.1f%%)", adj_var[2]), color = "Tissue", shape = "Tissue") +
  theme_portfolio() + theme(legend.position = "top")
save_figure(f06, "F06_visualization_only_adjusted_PCA", 9, 7)
write_tsv(adjusted_scores, "results/tables/figure_F06_source.tsv")
write_caption("F06_visualization_only_adjusted_PCA", "This transformed matrix is for visualization only. Differential expression uses the deposited matrix with run in the design.")

ann_df <- data.frame(Tissue = factor(metadata$tissue_group, levels = c("normal", "cancer")),
                     Run = factor(metadata$processing_run_proxy), row.names = metadata$gsm_id)
run_colors <- portfolio_run_colors(levels(ann_df$Run))
ann_colors <- list(Tissue = group_colors, Run = run_colors)
ph <- pheatmap(sample_correlation, color = portfolio_sequential(100),
               clustering_distance_rows = as.dist(1 - sample_correlation),
               clustering_distance_cols = as.dist(1 - sample_correlation),
               clustering_method = "average", annotation_row = ann_df, annotation_col = ann_df,
               annotation_colors = ann_colors, show_rownames = FALSE, show_colnames = FALSE,
               border_color = NA, main = "Sample Pearson correlation", fontsize = 8.5,
               fontsize_row = 7, fontsize_col = 7, silent = TRUE)
draw_f07 <- function() { grid::grid.newpage(); grid::grid.draw(ph$gtable) }
save_base_figure("F07_sample_correlation_heatmap", draw_f07, 10, 9)
write_tsv(as.data.table(sample_correlation, keep.rownames = "gsm_id"), "results/tables/figure_F07_source.tsv")
write_caption("F07_sample_correlation_heatmap", "Pearson correlations use all eligible representative genes; average-linkage clustering uses distance 1-r.")

message(sprintf("[02] PASS annotation: %.1f%% eligible probes; %d unique Entrez genes", 100 * eligible_fraction, nrow(gene_matrix)))
if (any(sample_metrics$candidate_exclusion_review)) {
  message(sprintf("[02] PAUSE/REVIEW: %d samples triggered at least two processed-data QC signals; no sample was excluded automatically.",
                  sum(sample_metrics$candidate_exclusion_review)))
} else {
  message("[02] QC: no sample triggered two independent processed-data QC signals")
}

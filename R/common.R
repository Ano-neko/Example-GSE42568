suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
  library(ggplot2)
  library(gridExtra)
  library(limma)
  library(AnnotationDbi)
})

options(stringsAsFactors = FALSE, width = 120)

parse_project_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match("--project", args)
  if (is.na(idx) || idx == length(args)) stop("Required argument: --project PATH")
  normalizePath(args[idx + 1], mustWork = TRUE)
}

init_project <- function() {
  project <- parse_project_arg()
  setwd(project)
  cfg <- yaml::read_yaml("config/analysis.yml")
  set.seed(as.integer(cfg$project$seed))
  Sys.setenv(TZ = cfg$project$timezone)
  list(project = project, cfg = cfg)
}

ensure_dirs <- function() {
  dirs <- c("data/metadata", "data/derived", "results/tables", "results/figures",
            "results/logs", "results/checksums", "reports", "docs", "tmp/pdfs", "output/pdf")
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
  invisible(path)
}

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub(" .*", "", out[[1]])
}

strip_geo_quotes <- function(x) {
  x <- sub('^"', '', x)
  sub('"$', '', x)
}

read_geo_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    out <- c(out, line)
    if (identical(line, "!series_matrix_table_begin")) break
  }
  out
}

geo_sample_rows <- function(header, key) {
  prefix <- paste0("!", key, "\t")
  lines <- header[startsWith(header, prefix)]
  lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
    strip_geo_quotes(fields)
  })
}

read_expression_matrix <- function(path, expected_rows) {
  cmd <- sprintf("gzip -dc %s", shQuote(path))
  dt <- data.table::fread(cmd = cmd, skip = "!series_matrix_table_begin",
                          nrows = expected_rows, check.names = FALSE,
                          data.table = TRUE, showProgress = FALSE)
  probe_ids <- dt[[1]]
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "double"
  rownames(mat) <- probe_ids
  mat
}

robust_z <- function(x) {
  center <- median(x, na.rm = TRUE)
  scale <- mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) return(rep(0, length(x)))
  (x - center) / scale
}

save_figure <- function(plot, stem, width = 10, height = 7, dpi = 600) {
  pdf_path <- file.path("results/figures", paste0(stem, ".pdf"))
  png_path <- file.path("results/figures", paste0(stem, ".png"))
  svg_path <- file.path("results/figures", paste0(stem, ".svg"))
  ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(svg_path, plot = plot, width = width, height = height, device = svglite::svglite, bg = "white")
  invisible(c(pdf_path, png_path, svg_path))
}

save_base_figure <- function(stem, draw_fun, width = 10, height = 7, dpi = 600) {
  pdf_path <- file.path("results/figures", paste0(stem, ".pdf"))
  png_path <- file.path("results/figures", paste0(stem, ".png"))
  svg_path <- file.path("results/figures", paste0(stem, ".svg"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height)
  draw_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  draw_fun()
  grDevices::dev.off()
  svglite::svglite(svg_path, width = width, height = height)
  draw_fun()
  grDevices::dev.off()
  invisible(c(pdf_path, png_path, svg_path))
}

write_caption <- function(stem, text) {
  writeLines(text, file.path("results/figures", paste0(stem, ".caption.txt")), useBytes = TRUE)
}

theme_portfolio <- function() {
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "#2F3A45"),
      axis.ticks = element_line(linewidth = 0.4, colour = "#2F3A45"),
      axis.ticks.length = grid::unit(2.2, "pt"),
      axis.title = element_text(size = 10, colour = "#2F3A45"),
      axis.text = element_text(size = 8.5, colour = "#55606C"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 13, colour = "#20262E", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 9.2, colour = "#667085", margin = margin(b = 8)),
      plot.margin = margin(9, 12, 9, 9),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9, colour = "#2F3A45"),
      legend.text = element_text(size = 8.2, colour = "#44505C"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9.5, colour = "#2F3A45"),
      panel.spacing = grid::unit(1.1, "lines")
    )
}

portfolio_colors <- c(
  cancer = "#FF8899",
  normal = "#7799CC",
  stable = "#779977",
  cancer_dark = "#D7657A",
  normal_dark = "#536F9F",
  stable_dark = "#5F7C5F",
  cancer_light = "#FFD6DD",
  normal_light = "#D8E3F2",
  stable_light = "#DCE7DC",
  neutral_dark = "#2F3A45",
  neutral_mid = "#7A8490",
  neutral_light = "#D9DEE5",
  neutral_pale = "#F4F6F8"
)

group_colors <- portfolio_colors[c("normal", "cancer")]
group_shapes <- c(normal = 17, cancer = 16)

portfolio_run_colors <- function(levels) {
  setNames(grDevices::hcl.colors(length(levels), "Set 3"), levels)
}

portfolio_diverging <- function(n = 101L) {
  grDevices::colorRampPalette(c(portfolio_colors[["normal_dark"]], "#FFFFFF", portfolio_colors[["cancer_dark"]]))(n)
}

portfolio_sequential <- function(n = 100L) {
  grDevices::colorRampPalette(c("#F5F7FA", portfolio_colors[["normal_light"]], portfolio_colors[["normal"]], portfolio_colors[["normal_dark"]]))(n)
}

fit_limma_model <- function(mat, metadata, include_run = TRUE, trend = TRUE, robust = TRUE,
                            treat_fc = 1.2) {
  metadata <- as.data.frame(metadata)
  metadata$tissue_group <- factor(metadata$tissue_group, levels = c("normal", "cancer"))
  if (include_run) {
    metadata$processing_run_proxy <- factor(metadata$processing_run_proxy)
    design <- model.matrix(~ processing_run_proxy + tissue_group, data = metadata)
  } else {
    design <- model.matrix(~ tissue_group, data = metadata)
  }
  if (qr(design)$rank != ncol(design)) stop("Design matrix is rank deficient")
  coef_name <- "tissue_groupcancer"
  if (!coef_name %in% colnames(design)) stop("Cancer-minus-normal coefficient is not estimable")
  fit0 <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit0, trend = trend, robust = robust)
  fit_treat <- limma::treat(fit0, lfc = log2(treat_fc), trend = trend, robust = robust)
  list(fit = fit, treat = fit_treat, design = design, coef_name = coef_name,
       coef_index = match(coef_name, colnames(design)))
}

extract_limma <- function(model) {
  j <- model$coef_index
  fit <- model$fit
  tr <- model$treat
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = fit$df.total)
  data.table(
    feature_id = rownames(fit$coefficients),
    log2FC = fit$coefficients[, j],
    fold_change_cancer_over_normal = 2 ^ fit$coefficients[, j],
    average_expression = fit$Amean,
    standard_error = se,
    ci95_low = fit$coefficients[, j] - crit * se,
    ci95_high = fit$coefficients[, j] + crit * se,
    moderated_t = fit$t[, j],
    p_value = fit$p.value[, j],
    bh_fdr = p.adjust(fit$p.value[, j], method = "BH"),
    treat_t = tr$t[, j],
    treat_p_value = tr$p.value[, j],
    treat_bh_fdr = p.adjust(tr$p.value[, j], method = "BH")
  )
}

hypergeom_ora <- function(selected, universe, term2gene, term_names, min_size, max_size) {
  selected <- intersect(unique(as.character(selected)), unique(as.character(universe)))
  term2gene <- unique(term2gene[gene_id %in% universe, .(term_id, gene_id)])
  sizes <- term2gene[, .(set_size = uniqueN(gene_id)), by = term_id]
  keep <- sizes[set_size >= min_size & set_size <= max_size, term_id]
  term2gene <- term2gene[term_id %in% keep]
  n <- length(selected)
  N <- length(unique(universe))
  ans <- term2gene[, {
    genes <- unique(gene_id)
    hit <- intersect(genes, selected)
    k <- length(hit); K <- length(genes)
    list(hit_count = k, selected_count = n, set_size = K, universe_size = N,
         p_value = if (n > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA_real_,
         hit_gene_ids = paste(sort(hit), collapse = ";"))
  }, by = term_id]
  ans[, term_name := unname(term_names[term_id])]
  ans[, bh_fdr := p.adjust(p_value, method = "BH")]
  ans[, gene_ratio := ifelse(selected_count > 0, hit_count / selected_count, NA_real_)]
  setorder(ans, p_value, -hit_count, term_id)
  ans
}

source("R/common.R")
ctx <- init_project()
cfg <- ctx$cfg
ensure_dirs()

message("[01] Validating official inputs and building sample metadata")
matrix_path <- cfg$input$matrix_file
platform_path <- cfg$input$platform_annotation_file
stopifnot(file.exists(matrix_path), file.exists(platform_path))

observed_matrix_sha <- sha256_file(matrix_path)
if (!identical(observed_matrix_sha, cfg$input$expected_matrix_sha256)) {
  stop(sprintf("Matrix SHA-256 mismatch: expected %s, observed %s",
               cfg$input$expected_matrix_sha256, observed_matrix_sha))
}

gzip_ok <- function(path) {
  status <- system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE)
  identical(status, 0L)
}
if (!gzip_ok(matrix_path) || !gzip_ok(platform_path)) stop("A gzip integrity test failed")

header <- read_geo_header(matrix_path)
get_one <- function(key) {
  x <- geo_sample_rows(header, key)
  if (!length(x)) return(NULL)
  x[[1]]
}

gsm <- get_one("Sample_geo_accession")
titles <- get_one("Sample_title")
sources <- get_one("Sample_source_name_ch1")
platforms <- get_one("Sample_platform_id")
row_counts <- suppressWarnings(as.integer(get_one("Sample_data_row_count")))
sample_status <- get_one("Sample_status")
sample_type <- get_one("Sample_type")
data_processing <- get_one("Sample_data_processing")

n_expected <- as.integer(cfg$input$expected_samples)
required_vectors <- list(gsm = gsm, title = titles, source = sources, platform = platforms)
bad_lengths <- names(required_vectors)[vapply(required_vectors, length, integer(1)) != n_expected]
if (length(bad_lengths)) stop("Unexpected GEO header vector lengths: ", paste(bad_lengths, collapse = ", "))

characteristic_rows <- geo_sample_rows(header, "Sample_characteristics_ch1")
if (!length(characteristic_rows)) stop("No sample characteristics found in series matrix header")
if (any(vapply(characteristic_rows, length, integer(1)) != n_expected)) {
  stop("At least one sample-characteristics row does not have 121 values")
}

normalise_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

chars <- list()
for (row in characteristic_rows) {
  pieces <- strsplit(row, ":", fixed = TRUE)
  keys_row <- vapply(pieces, function(z) normalise_key(z[[1]]), character(1))
  if (length(unique(keys_row)) != 1L) stop("A characteristics row contains inconsistent field names")
  key <- keys_row[[1]]
  values <- vapply(pieces, function(z) {
    if (length(z) < 2) return(NA_character_)
    val <- trimws(paste(z[-1], collapse = ":"))
    if (toupper(val) == "NA" || val == "") NA_character_ else val
  }, character(1))
  chars[[key]] <- values
}

metadata <- data.table(
  gsm_id = gsm,
  original_title = titles,
  original_source_name = sources,
  platform_id = platforms,
  row_count = row_counts,
  sample_status = if (length(sample_status)) sample_status else NA_character_,
  sample_type = if (length(sample_type)) sample_type else NA_character_,
  tissue_group = fifelse(grepl("normal", sources, ignore.case = TRUE), "normal",
                         fifelse(grepl("cancer", sources, ignore.case = TRUE), "cancer", NA_character_))
)

for (nm in names(chars)) metadata[[nm]] <- chars[[nm]]

numeric_candidates <- c("age", "size", "grade", "er_status", "lymph_node_status",
                        "relapse_free_survival_event", "overall_survival_event",
                        "relapse_free_survival_months", "overall_survival_months")
for (nm in intersect(numeric_candidates, names(metadata))) {
  metadata[[nm]] <- suppressWarnings(as.numeric(metadata[[nm]]))
}

parse_title_date <- function(title) {
  match <- regexec("_([0-9]{1,2})_([0-9]{1,2})_([0-9]{2})$", title)
  hit <- regmatches(title, match)[[1]]
  if (length(hit) != 4L) return(NA_character_)
  sprintf("20%02d-%02d-%02d", as.integer(hit[[4]]), as.integer(hit[[3]]), as.integer(hit[[2]]))
}
metadata[, processing_run_proxy := vapply(original_title, parse_title_date, character(1))]
metadata[, processing_run_source := fifelse(is.na(processing_run_proxy),
                                             "manual_raw_CEL_header_audit_from_specification",
                                             "sample_title_date_suffix")]
manual_map <- unlist(cfg$design$manual_run_proxy, use.names = TRUE)
for (id in names(manual_map)) metadata[gsm_id == id, processing_run_proxy := as.character(manual_map[[id]])]

metadata[, paired_identifier := NA_character_]
metadata[, normal_provenance := fifelse(tissue_group == "normal", "not_reported", "not_applicable")]
metadata[, included_primary := TRUE]
metadata[, exclusion_reason := NA_character_]

if (anyNA(metadata$tissue_group)) stop("Unresolved tissue group")
if (anyNA(metadata$processing_run_proxy)) stop("Unresolved processing-run proxy")
if (anyDuplicated(metadata$gsm_id)) stop("Duplicate GSM IDs in metadata")
if (metadata[, sum(tissue_group == "cancer")] != cfg$input$expected_cancer) stop("Cancer count mismatch")
if (metadata[, sum(tissue_group == "normal")] != cfg$input$expected_normal) stop("Normal count mismatch")
if (!all(metadata$platform_id == cfg$project$platform)) stop("Platform mismatch")

expr <- read_expression_matrix(matrix_path, as.integer(cfg$input$expected_probes))
if (!identical(colnames(expr), metadata$gsm_id)) stop("Matrix columns do not match GEO sample order")
if (nrow(expr) != cfg$input$expected_probes || ncol(expr) != cfg$input$expected_samples) stop("Matrix dimension mismatch")
if (anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) stop("Duplicate matrix IDs")
if (any(!is.finite(expr))) stop("Matrix contains missing or non-finite values")
if (max(expr) > 50 || min(expr) < -10) stop("Matrix scale is inconsistent with deposited log2 GC-RMA values")

constant <- apply(expr, 1, function(x) max(x) == min(x))
integrity <- data.table(
  metric = c("matrix_sha256", "rows_probe_sets", "columns_samples", "unique_probe_ids",
             "unique_sample_ids", "nonfinite_values", "minimum", "maximum",
             "constant_probe_sets", "AFFX_control_probe_sets", "cancer_samples",
             "normal_samples", "platform", "scale_assessment", "gzip_integrity"),
  value = c(observed_matrix_sha, nrow(expr), ncol(expr), uniqueN(rownames(expr)),
            uniqueN(colnames(expr)), sum(!is.finite(expr)), sprintf("%.10f", min(expr)),
            sprintf("%.10f", max(expr)), sum(constant), sum(startsWith(rownames(expr), "AFFX-")),
            metadata[, sum(tissue_group == "cancer")], metadata[, sum(tissue_group == "normal")],
            unique(metadata$platform_id), "consistent_with_log2_GCRMA", "passed")
)

run_counts <- dcast(metadata, processing_run_proxy ~ tissue_group, fun.aggregate = length, value.var = "gsm_id")
if (!"normal" %in% names(run_counts)) run_counts[, normal := 0L]
if (!"cancer" %in% names(run_counts)) run_counts[, cancer := 0L]
run_counts[, mixed_run := normal > 0 & cancer > 0]
setorder(run_counts, processing_run_proxy)

model_metadata <- copy(metadata)
model_metadata[, tissue_group := factor(tissue_group, levels = c("normal", "cancer"))]
model_metadata[, processing_run_proxy := factor(processing_run_proxy)]
design <- model.matrix(~ processing_run_proxy + tissue_group, data = model_metadata)
design_summary <- data.table(
  metric = c("design_formula", "design_rows", "design_columns", "design_rank",
             "full_rank", "residual_degrees_of_freedom", "run_levels",
             "mixed_run_levels", "mixed_run_samples", "contrast", "contrast_estimable"),
  value = c("expression ~ processing_run_proxy + tissue_group", nrow(design), ncol(design),
            qr(design)$rank, qr(design)$rank == ncol(design), nrow(design) - qr(design)$rank,
            uniqueN(metadata$processing_run_proxy), sum(run_counts$mixed_run),
            metadata[processing_run_proxy %in% run_counts[mixed_run == TRUE, processing_run_proxy], .N],
            "cancer-minus-normal", "tissue_groupcancer" %in% colnames(design))
)
if (qr(design)$rank != ncol(design) || !"tissue_groupcancer" %in% colnames(design)) {
  stop("Hard stop: cancer-minus-normal is not estimable after run adjustment")
}

age_lt50 <- metadata[tissue_group == "cancer" & age < 50, .N]
age_ge50 <- metadata[tissue_group == "cancer" & age >= 50, .N]
issues <- data.table(
  issue_id = c("META-001", "META-002", "META-003", "META-004", "ACQ-001"),
  severity = c("document", "limitation", "limitation", "document", "document"),
  issue = c(
    sprintf("Series narrative reports 20 tumors under 50, but sample-level metadata contain %d under 50 and %d at least 50.", age_lt50, age_ge50),
    "The provenance of the 17 normal breast tissues is not reported sufficiently to call them healthy, adjacent, or matched controls.",
    "Normal-sample age and pairing identifiers are absent; age adjustment and paired analysis are not supportable.",
    "Three tumor titles lack a date suffix; the frozen specification assigns 2004-11-26 after a limited raw-CEL header audit.",
    "Family SOFT was not retained in v1 because all required sample fields are present in the locked series-matrix header; this is a documented infrastructure omission, not a scientific-data substitution."
  ),
  resolution = c("Use sample-level values and disclose discrepancy", "Use the term normal breast tissue only",
                 "Use an unpaired design without age adjustment", "Retain manual map and provenance field",
                 "Reproducibility relies on the locked series matrix plus official accession")
)

field_dictionary <- data.table(
  field = names(metadata),
  type = vapply(metadata, function(x) class(x)[1], character(1)),
  role = fifelse(names(metadata) == "tissue_group", "primary_exposure",
          fifelse(names(metadata) == "processing_run_proxy", "fixed_effect_batch_proxy",
          fifelse(names(metadata) %in% c("gsm_id", "original_title", "platform_id"), "identity_or_provenance",
          fifelse(names(metadata) %in% numeric_candidates, "tumor_only_descriptive_not_primary_covariate", "metadata")))),
  missing_value_policy = "NA is explicit; not_reported and not_applicable are retained in companion fields or documentation"
)

header_value <- function(path, pattern) {
  if (!file.exists(path)) return(NA_character_)
  lines <- readLines(path, warn = FALSE)
  hit <- grep(pattern, lines, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) NA_character_ else trimws(sub("^[^:]+:", "", tail(hit, 1)))
}
source_manifest <- data.table(
  accession = c("GSE42568", "GPL570", "GSE42568", "Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor"),
  asset_role = c("processed_series_matrix", "platform_annotation_snapshot", "family_SOFT_optional_not_retained",
                 "probe_annotation_package", "human_gene_annotation_package", "GO_annotation_package", "Reactome_annotation_package"),
  source_url = c(cfg$input$matrix_url, cfg$input$platform_annotation_url,
                 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/soft/GSE42568_family.soft.gz",
                 "https://bioconductor.org/packages/hgu133plus2.db",
                 "https://bioconductor.org/packages/org.Hs.eg.db",
                 "https://bioconductor.org/packages/GO.db",
                 "https://bioconductor.org/packages/reactome.db"),
  local_path = c(matrix_path, platform_path, NA_character_, "R_library:hgu133plus2.db", "R_library:org.Hs.eg.db",
                 "R_library:GO.db", "R_library:reactome.db"),
  retrieved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  http_last_modified = c(header_value("data/metadata/GSE42568_series_matrix.headers.txt", "Last-Modified"),
                         header_value("data/metadata/GPL570_annotation.headers.txt", "Last-Modified"), rep(NA_character_, 5)),
  content_length_bytes = c(file.info(matrix_path)$size, file.info(platform_path)$size, NA_real_, rep(NA_real_, 4)),
  sha256 = c(observed_matrix_sha, sha256_file(platform_path), NA_character_, rep(NA_character_, 4)),
  compression_test = c("passed", "passed", "not_retained", rep("not_applicable", 4)),
  parser_or_package_version = c(as.character(packageVersion("data.table")), "not_used_for_primary_mapping", "not_applicable",
                                as.character(packageVersion("hgu133plus2.db")), as.character(packageVersion("org.Hs.eg.db")),
                                as.character(packageVersion("GO.db")), as.character(packageVersion("reactome.db"))),
  license_or_terms_url = c("https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html",
                           "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html",
                           "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html",
                           rep("https://bioconductor.org/about/licenses/", 4)),
  status = c("locked", "locked_auxiliary", "documented_omission", rep("installed_and_locked", 4))
)

write_tsv(metadata, "data/metadata/sample_manifest.tsv")
write_tsv(field_dictionary, "data/metadata/field_dictionary.tsv")
write_tsv(issues, "data/metadata/metadata_issues.tsv")
write_tsv(source_manifest, "data/metadata/source_manifest.tsv")
write_tsv(integrity, "results/tables/matrix_integrity.tsv")
write_tsv(run_counts, "results/tables/run_group_counts.tsv")
write_tsv(design_summary, "results/tables/design_diagnostics.tsv")
saveRDS(expr, "data/derived/expression_probe_matrix.rds", compress = "xz")
saveRDS(metadata, "data/derived/sample_metadata.rds", compress = "xz")
saveRDS(header, "data/derived/geo_matrix_header.rds", compress = "xz")

message(sprintf("[01] PASS: %d probes x %d samples; design rank %d/%d; %d mixed-run samples",
                nrow(expr), ncol(expr), qr(design)$rank, ncol(design),
                metadata[processing_run_proxy %in% run_counts[mixed_run == TRUE, processing_run_proxy], .N]))

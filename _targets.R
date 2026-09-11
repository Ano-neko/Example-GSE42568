library(targets)

project <- normalizePath(".", mustWork = TRUE)
resolve_runtime <- function(command, override_var) {
  override <- Sys.getenv(override_var, unset = "")
  prefix <- Sys.getenv("BIOINFO_ENV_PREFIX", unset = "")
  candidate <- if (nzchar(override)) {
    override
  } else if (nzchar(prefix)) {
    file.path(prefix, "bin", command)
  } else {
    unname(Sys.which(command))
  }
  if (nzchar(candidate) && !grepl("[/\\\\]", candidate)) {
    candidate <- unname(Sys.which(candidate))
  }
  if (!nzchar(candidate) || !file.exists(candidate)) {
    stop("Cannot resolve ", command, "; activate the Conda environment or set BIOINFO_ENV_PREFIX")
  }
  normalizePath(candidate, mustWork = TRUE)
}
rscript <- resolve_runtime("Rscript", "R_BIN")
python <- resolve_runtime("python", "PY_BIN")
run_r <- function(script) {
  status <- system2(rscript, c(script, "--project", project))
  if (!identical(status, 0L)) stop(script, " failed with status ", status)
}

list(
  tar_target(raw_inputs, {
    status <- system2("bash", c("scripts/download_inputs.sh", project))
    if (!identical(status, 0L)) stop("download failed")
    c("data/raw/GSE42568_series_matrix.txt.gz", "data/raw/GPL570.annot.gz")
  }, format = "file"),
  tar_target(stage01, { raw_inputs; run_r("R/01_acquire_validate.R"); "data/derived/expression_probe_matrix.rds" }, format = "file"),
  tar_target(stage02, { stage01; run_r("R/02_qc_annotation.R"); "data/derived/expression_gene_representative.rds" }, format = "file"),
  tar_target(qc_review, { stage02; run_r("R/02b_qc_review.R"); "results/tables/sample_exclusion_ledger.tsv" }, format = "file"),
  tar_target(stage03, { qc_review; run_r("R/03_de_sensitivity.R"); "results/tables/de_primary_all_genes.tsv" }, format = "file"),
  tar_target(stage04, { stage03; run_r("R/04_enrichment.R"); "results/tables/gsea_all_gene_sets.tsv" }, format = "file"),
  tar_target(guards, { stage04; run_r("tests/run_tests.R"); "results/tables/reproducibility_guard_tests.tsv" }, format = "file"),
  tar_target(scientific_qa, { guards; run_r("R/05_scientific_qa.R"); "results/tables/acceptance_checks.tsv" }, format = "file"),
  tar_target(reports, {
    scientific_qa
    status <- system2(python, c("scripts/build_reports.py", "--project", project))
    if (!identical(status, 0L)) stop("report build failed")
    c("output/pdf/client_report.pdf", "output/pdf/methods_appendix.pdf", "output/pdf/qc_appendix.pdf")
  }, format = "file"),
  tar_target(public_release_qa, {
    reports
    status <- system2(python, c("scripts/validate_public_release.py", "--project", project))
    if (!identical(status, 0L)) stop("public release QA failed")
    "results/tables/public_portfolio_acceptance.tsv"
  }, format = "file")
)

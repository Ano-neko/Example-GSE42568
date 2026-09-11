source("R/common.R")
ctx <- init_project()
ensure_dirs()

message("[test] Running reproducibility and scientific-guard fixtures")
tests <- list()
record <- function(id, condition, detail) {
  tests[[length(tests) + 1L]] <<- data.table(test_id = id, passed = isTRUE(condition), detail = detail)
}

# A tiny balanced fixture locks the interpretation of the coefficient sign.
fixture_md <- data.frame(
  tissue_group = c("normal", "cancer", "normal", "cancer"),
  processing_run_proxy = c("r1", "r1", "r2", "r2")
)
fixture_y <- matrix(c(5, 7, 6, 8, 9, 9, 9, 9), nrow = 2, byrow = TRUE,
                    dimnames = list(c("positive_control", "null_control"), paste0("s", 1:4)))
fixture_fit <- fit_limma_model(fixture_y, fixture_md, include_run = TRUE, trend = FALSE, robust = FALSE)
coef_positive <- fixture_fit$fit$coefficients["positive_control", fixture_fit$coef_index]
coef_null <- fixture_fit$fit$coefficients["null_control", fixture_fit$coef_index]
record("T01_contrast_cancer_minus_normal_positive", abs(coef_positive - 2) < 1e-10,
       sprintf("Expected +2; observed %.12f", coef_positive))
record("T02_contrast_null_control", abs(coef_null) < 1e-10,
       sprintf("Expected 0; observed %.12f", coef_null))

manifest <- fread("data/metadata/sample_manifest.tsv")
record("T03_sample_groups_locked", nrow(manifest) == 121L && sum(manifest$tissue_group == "cancer") == 104L &&
         sum(manifest$tissue_group == "normal") == 17L,
       "Expected 121 samples: 104 cancer and 17 normal")
paired_nonmissing <- if ("paired_identifier" %in% names(manifest)) sum(!is.na(manifest$paired_identifier) & manifest$paired_identifier != "") else NA_integer_
record("T04_unpaired_design_guard", identical(paired_nonmissing, 0L),
       sprintf("Non-missing pairing identifiers: %s", paired_nonmissing))

probe_audit <- fread("results/tables/probe_mapping_audit.tsv")
record("T05_one_to_many_not_primary", all(probe_audit[mapping_status == "one_to_many", primary_eligible == FALSE]),
       "Every one-to-many probe is excluded from primary annotation")
record("T06_x_at_not_primary", all(probe_audit[suffix_class == "x_at", primary_eligible == FALSE]),
       "Every _x_at probe is excluded from primary annotation")
record("T07_controls_not_primary", all(probe_audit[is_affx_control == TRUE, primary_eligible == FALSE]),
       "Every AFFX control is excluded from primary annotation")

selection <- fread("results/tables/probe_selection_audit.tsv", colClasses = list(character = "entrez_id"))
selected <- selection[selected_representative == TRUE]
record("T08_one_representative_per_gene", nrow(selected) == uniqueN(selection$entrez_id) &&
         all(selected[, .N, by = entrez_id]$N == 1L),
       sprintf("Selected %d representatives for %d Entrez genes", nrow(selected), uniqueN(selection$entrez_id)))

de <- fread("results/tables/de_primary_all_genes.tsv", colClasses = list(character = "entrez_id"))
universe <- fread("results/tables/enrichment_universe_membership.tsv", colClasses = list(character = "entrez_id"))
record("T09_universes_subset_tested", all(universe$entrez_id %in% de$entrez_id),
       "Every resource-universe member occurs in the complete P1 tested-gene table")

gsea_rank <- fread("results/tables/gsea_rank_vector.tsv", colClasses = list(character = "entrez_id"))
record("T10_gsea_rank_complete_unique", nrow(gsea_rank) == nrow(de) && !anyDuplicated(gsea_rank$entrez_id) &&
         setequal(gsea_rank$entrez_id, de$entrez_id),
       "GSEA rank contains every tested gene exactly once")
record("T11_gsea_rank_descending", !is.unsorted(-gsea_rank$moderated_t),
       "Moderated t is monotonically descending; exact ties use Entrez ID")

license_manifest <- fread("results/tables/gene_set_manifest.tsv")
record("T12_kegg_license_gate", license_manifest[collection == "KEGG_ORA", status] == "not_run_license_gate",
       "KEGG was not run because commercial permission was not documented")
record("T13_msigdb_license_gate", license_manifest[collection == "MSigDB_Hallmark_GSEA", status] == "not_run_license_review_unconfirmed",
       "MSigDB Hallmark was not run because commercial-use review was unconfirmed")

test_dt <- rbindlist(tests)
write_tsv(test_dt, "results/tables/reproducibility_guard_tests.tsv")
if (any(!test_dt$passed)) {
  print(test_dt[passed == FALSE])
  stop("One or more reproducibility/scientific guard tests failed")
}
message(sprintf("[test] PASS: %d/%d guards", sum(test_dt$passed), nrow(test_dt)))

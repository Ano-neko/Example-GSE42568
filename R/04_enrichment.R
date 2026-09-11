source("R/common.R")
ctx <- init_project()
cfg <- ctx$cfg
ensure_dirs()

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(reactome.db)
  library(fgsea)
  library(gridExtra)
})

message("[04] Licensed enrichment with tested-and-annotated universes")
de <- fread("results/tables/de_primary_all_genes.tsv", colClasses = list(character = "entrez_id"))
licenses <- yaml::read_yaml("config/gene_set_licenses.yml")
stopifnot(nrow(de) == 20357L, !anyDuplicated(de$entrez_id))
stopifnot(identical(licenses$KEGG$status, "not_run_license_gate"))
stopifnot(identical(licenses$MSigDB_Hallmark$status, "not_run_license_review_unconfirmed"))

ora_min <- as.integer(cfg$enrichment$ora_min_size)
ora_max <- as.integer(cfg$enrichment$ora_max_size)
gsea_min <- as.integer(cfg$enrichment$gsea_min_size)
gsea_max <- as.integer(cfg$enrichment$gsea_max_size)
gsea_exp <- as.numeric(cfg$enrichment$gsea_score_exponent)
fdr_cut <- as.numeric(cfg$enrichment$bh_fdr)
tested <- de$entrez_id
priority_up <- de[priority_effect == TRUE & log2FC > 0, entrez_id]
priority_down <- de[priority_effect == TRUE & log2FC < 0, entrez_id]
symbol_map <- setNames(de$symbol_at_run, de$entrez_id)

add_hit_symbols <- function(x) {
  x[, hit_gene_symbols := vapply(strsplit(hit_gene_ids, ";", fixed = TRUE), function(ids) {
    if (length(ids) == 1L && ids == "") return("")
    paste(unname(symbol_map[ids]), collapse = ";")
  }, character(1))]
  x
}

# GO mapping uses propagated GOALL annotations, while universes remain restricted
# to genes actually tested by P1 and annotated in each ontology.
message("[04] Building GOALL mapping")
go_sel <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = tested, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")
))
go_map <- unique(as.data.table(go_sel)[!is.na(GOALL) & ONTOLOGYALL %in% c("BP", "CC", "MF"),
                                     .(term_id = GOALL, gene_id = ENTREZID, ontology = ONTOLOGYALL)])
go_terms <- unique(go_map$term_id)
go_names_dt <- as.data.table(AnnotationDbi::select(GO.db, keys = go_terms, keytype = "GOID", columns = "TERM"))
go_names_dt <- unique(go_names_dt[!is.na(TERM), .(term_id = GOID, term_name = TERM)])
go_name_vec <- setNames(go_names_dt$term_name, go_names_dt$term_id)

# Original GO set sizes are package-wide; no untested gene is added to the ORA universe.
go_all_list <- as.list(org.Hs.egGO2ALLEGS)
go_original_size <- vapply(go_all_list[go_terms], function(x) length(unique(unname(x))), integer(1))

go_results <- list()
go_universe_rows <- list()
go_manifest_rows <- list()
counter <- 1L
for (ont in c("BP", "CC", "MF")) {
  t2g <- go_map[ontology == ont, .(term_id, gene_id)]
  universe <- sort(unique(t2g$gene_id))
  go_universe_rows[[ont]] <- data.table(resource = paste0("GO_", ont), entrez_id = universe)
  sizes <- t2g[, .(intersected_size = uniqueN(gene_id)), by = term_id]
  eligible_sets <- sizes[intersected_size >= ora_min & intersected_size <= ora_max, .N]
  go_manifest_rows[[ont]] <- data.table(
    collection = paste0("GO_", ont, "_ORA"), status = "completed", role = ifelse(ont == "BP", "primary_ORA", "supplement_ORA"),
    source = "org.Hs.eg.db GOALL plus GO.db", source_version = paste0("org.Hs.eg.db ", packageVersion("org.Hs.eg.db"), "; GO.db ", packageVersion("GO.db")),
    license_status = licenses$GO$status, universe_definition = "P1-tested Entrez genes with at least one annotation in this ontology",
    universe_size = length(universe), raw_gene_sets = uniqueN(t2g$term_id), eligible_gene_sets = eligible_sets,
    min_size = ora_min, max_size = ora_max, ranking = NA_character_
  )
  for (direction in c("cancer_up", "normal_up")) {
    selected0 <- if (direction == "cancer_up") priority_up else priority_down
    selected <- intersect(selected0, universe)
    if (length(selected) < 10L) {
      go_results[[counter]] <- data.table(
        status = "not_evaluable_fewer_than_10_mapped_priority_genes", resource = "GO", ontology = ont,
        direction = direction, term_id = NA_character_, term_name = NA_character_, selected_count = length(selected),
        universe_size = length(universe), p_value = NA_real_, bh_fdr = NA_real_
      )
    } else {
      ans <- hypergeom_ora(selected, universe, t2g, go_name_vec, ora_min, ora_max)
      ans[, `:=`(status = "completed", resource = "GO", ontology = ont, direction = direction,
                 original_set_size = unname(go_original_size[term_id]), intersected_set_size = set_size,
                 selected_gene_ids = paste(sort(selected), collapse = ";"),
                 annotation_snapshot = paste0("org.Hs.eg.db_", packageVersion("org.Hs.eg.db"), "__GO.db_", packageVersion("GO.db")))]
      ans <- add_hit_symbols(ans)
      go_results[[counter]] <- ans
    }
    counter <- counter + 1L
  }
}
go_ora <- rbindlist(go_results, fill = TRUE)
go_ora <- go_ora[order(ontology, direction, bh_fdr, p_value, term_id, na.last = TRUE)]
write_tsv(go_ora, "results/tables/go_ora_all_terms.tsv")
fwrite(
  go_ora,
  "results/tables/go_ora_all_terms.tsv.gz",
  sep = "\t",
  na = "NA",
  quote = FALSE,
  compress = "gzip"
)

# Reactome CC0 pathway mapping, restricted to human R-HSA pathways.
message("[04] Building Reactome mapping")
rx_gene_list <- as.list(reactomeEXTID2PATHID)
rx_gene_list <- rx_gene_list[intersect(names(rx_gene_list), tested)]
rx_map <- rbindlist(lapply(names(rx_gene_list), function(id) {
  paths <- unique(rx_gene_list[[id]])
  paths <- paths[grepl("^R-HSA-", paths)]
  if (!length(paths)) return(NULL)
  data.table(term_id = paths, gene_id = id)
}))
rx_map <- unique(rx_map)
rx_path_ids <- unique(rx_map$term_id)
rx_name_list <- as.list(reactomePATHID2NAME)
rx_name_vec <- vapply(rx_path_ids, function(id) {
  x <- rx_name_list[[id]]
  if (is.null(x) || !length(x)) id else sub("^Homo sapiens: ", "", x[[1]])
}, character(1))
names(rx_name_vec) <- rx_path_ids
rx_original_list <- as.list(reactomePATHID2EXTID)
rx_original_size <- vapply(rx_path_ids, function(id) length(unique(unname(rx_original_list[[id]]))), integer(1))
names(rx_original_size) <- rx_path_ids
rx_universe <- sort(unique(rx_map$gene_id))
rx_sizes <- rx_map[, .(intersected_size = uniqueN(gene_id)), by = term_id]

rx_ora_list <- lapply(c("cancer_up", "normal_up"), function(direction) {
  selected0 <- if (direction == "cancer_up") priority_up else priority_down
  selected <- intersect(selected0, rx_universe)
  if (length(selected) < 10L) return(data.table(
    status = "not_evaluable_fewer_than_10_mapped_priority_genes", resource = "Reactome", direction = direction,
    term_id = NA_character_, term_name = NA_character_, selected_count = length(selected), universe_size = length(rx_universe),
    p_value = NA_real_, bh_fdr = NA_real_))
  ans <- hypergeom_ora(selected, rx_universe, rx_map, rx_name_vec, ora_min, ora_max)
  ans[, `:=`(status = "completed", resource = "Reactome", direction = direction,
             original_set_size = unname(rx_original_size[term_id]), intersected_set_size = set_size,
             selected_gene_ids = paste(sort(selected), collapse = ";"),
             annotation_snapshot = paste0("reactome.db_", packageVersion("reactome.db")))]
  add_hit_symbols(ans)
})
rx_ora <- rbindlist(rx_ora_list, fill = TRUE)
kegg_gate <- data.table(
  status = "KEGG_NOT_RUN_LICENSE_GATE", resource = "KEGG", direction = "not_applicable",
  term_id = NA_character_, term_name = NA_character_, selected_count = NA_integer_, universe_size = NA_integer_,
  p_value = NA_real_, bh_fdr = NA_real_, annotation_snapshot = "No KEGG API, cache, or pathway definition accessed"
)
pathway_ora <- rbindlist(list(rx_ora, kegg_gate), fill = TRUE)
pathway_ora <- pathway_ora[order(resource, direction, bh_fdr, p_value, term_id, na.last = TRUE)]
write_tsv(pathway_ora, "results/tables/pathway_ora_all_terms.tsv")

# Reactome preranked GSEA fallback. Stable tie-breaking changes ordering only for
# equal t statistics; it does not alter the statistics.
rank_dt <- de[, .(entrez_id, moderated_t)]
setorder(rank_dt, -moderated_t, entrez_id)
tie_member <- duplicated(rank_dt$moderated_t) | duplicated(rank_dt$moderated_t, fromLast = TRUE)
tie_fraction <- mean(tie_member)
tie_unique_fraction <- 1 - uniqueN(rank_dt$moderated_t) / nrow(rank_dt)
write_tsv(rank_dt[, rank := .I], "results/tables/gsea_rank_vector.tsv")
if (tie_fraction > 0.01) stop("GSEA_PAUSE: more than 1% of tested genes participate in exact moderated-t ties")

rx_pathways <- split(rx_map$gene_id, rx_map$term_id)
ranks <- setNames(rank_dt$moderated_t, rank_dt$entrez_id)
set.seed(as.integer(cfg$project$seed))
worker_n <- min(8L, parallel::detectCores())
bp_param <- BiocParallel::MulticoreParam(workers = worker_n, progressbar = FALSE, RNGseed = as.integer(cfg$project$seed))
message(sprintf("[04] Running Reactome fgseaMultilevel on %d pathways with %d workers", length(rx_pathways), worker_n))
gsea <- suppressWarnings(fgseaMultilevel(
  pathways = rx_pathways, stats = ranks, minSize = gsea_min, maxSize = gsea_max,
  eps = 1e-50, scoreType = "std", gseaParam = gsea_exp, BPPARAM = bp_param
))
gsea_dt <- as.data.table(gsea)
setnames(gsea_dt, c("pathway", "pval", "padj", "size", "leadingEdge"),
         c("term_id", "p_value", "bh_fdr", "intersected_set_size", "leading_edge_ids"))
gsea_dt[, leading_edge_ids := vapply(leading_edge_ids, function(x) paste(x, collapse = ";"), character(1))]
gsea_dt[, `:=`(
  status = "completed", collection = "Reactome_preranked_GSEA_fallback", term_name = unname(rx_name_vec[term_id]),
  original_set_size = unname(rx_original_size[term_id]), direction = ifelse(NES > 0, "cancer_up", "normal_up"),
  rank_statistic = "P1 cancer-minus-normal moderated t", rank_gene_count = length(ranks),
  rank_tie_member_fraction = tie_fraction, score_exponent = gsea_exp,
  source_snapshot = paste0("reactome.db_", packageVersion("reactome.db"), "__fgsea_", packageVersion("fgsea"))
)]
gsea_dt[, leading_edge_symbols := vapply(strsplit(leading_edge_ids, ";", fixed = TRUE), function(ids) {
  if (length(ids) == 1L && ids == "") return("")
  paste(unname(symbol_map[ids]), collapse = ";")
}, character(1))]
msigdb_gate <- data.table(
  status = "MSIGDB_HALLMARK_NOT_RUN_LICENSE_REVIEW_UNCONFIRMED", collection = "MSigDB_Hallmark",
  term_id = NA_character_, term_name = NA_character_, p_value = NA_real_, bh_fdr = NA_real_, NES = NA_real_,
  direction = "not_applicable", rank_statistic = "not_run", rank_gene_count = length(ranks),
  source_snapshot = "No MSigDB collection downloaded or redistributed"
)
gsea_all <- rbindlist(list(gsea_dt, msigdb_gate), fill = TRUE)
gsea_all <- gsea_all[order(collection, bh_fdr, -abs(NES), term_id, na.last = TRUE)]
write_tsv(gsea_all, "results/tables/gsea_all_gene_sets.tsv")

universe_table <- rbindlist(c(go_universe_rows, list(Reactome = data.table(resource = "Reactome", entrez_id = rx_universe))))
write_tsv(universe_table, "results/tables/enrichment_universe_membership.tsv")

rx_manifest <- data.table(
  collection = c("Reactome_ORA", "Reactome_preranked_GSEA_fallback", "KEGG_ORA", "MSigDB_Hallmark_GSEA"),
  status = c("completed", "completed", "not_run_license_gate", "not_run_license_review_unconfirmed"),
  role = c("pathway_ORA_replacement", "primary_GSEA_fallback", "not_run", "not_run"),
  source = c("reactome.db", "reactome.db", "No KEGG resource accessed", "No MSigDB resource accessed"),
  source_version = c(as.character(packageVersion("reactome.db")), as.character(packageVersion("reactome.db")), NA_character_, NA_character_),
  license_status = c(licenses$Reactome$status, licenses$Reactome$status, licenses$KEGG$status, licenses$MSigDB_Hallmark$status),
  universe_definition = c(rep("P1-tested Entrez genes with at least one Reactome human pathway annotation", 2), NA_character_, NA_character_),
  universe_size = c(length(rx_universe), length(rx_universe), NA_integer_, NA_integer_),
  raw_gene_sets = c(uniqueN(rx_map$term_id), uniqueN(rx_map$term_id), NA_integer_, NA_integer_),
  eligible_gene_sets = c(rx_sizes[intersected_size >= ora_min & intersected_size <= ora_max, .N],
                         rx_sizes[intersected_size >= gsea_min & intersected_size <= gsea_max, .N], NA_integer_, NA_integer_),
  min_size = c(ora_min, gsea_min, NA_integer_, NA_integer_), max_size = c(ora_max, gsea_max, NA_integer_, NA_integer_),
  ranking = c(NA_character_, "P1 moderated t, descending; exact ties by Entrez ID", NA_character_, NA_character_)
)
gene_set_manifest <- rbindlist(c(go_manifest_rows, list(rx_manifest)), fill = TRUE)
write_tsv(gene_set_manifest, "results/tables/gene_set_manifest.tsv")

enrichment_diag <- data.table(
  metric = c("priority_cancer_up", "priority_normal_up", "go_bp_universe", "go_cc_universe", "go_mf_universe",
             "reactome_universe", "reactome_eligible_ora_sets", "reactome_eligible_gsea_sets",
             "gsea_rank_genes", "gsea_tie_member_fraction", "gsea_duplicate_statistic_fraction",
             "go_bp_significant_rows", "reactome_ora_significant_rows", "reactome_gsea_significant_rows",
             "kegg_status", "msigdb_hallmark_status"),
  value = c(length(priority_up), length(priority_down),
            nrow(go_universe_rows$BP), nrow(go_universe_rows$CC), nrow(go_universe_rows$MF), length(rx_universe),
            rx_sizes[intersected_size >= ora_min & intersected_size <= ora_max, .N],
            rx_sizes[intersected_size >= gsea_min & intersected_size <= gsea_max, .N], length(ranks),
            tie_fraction, tie_unique_fraction,
            go_ora[ontology == "BP" & bh_fdr < fdr_cut, .N], rx_ora[bh_fdr < fdr_cut, .N], gsea_dt[bh_fdr < fdr_cut, .N],
            licenses$KEGG$status, licenses$MSigDB_Hallmark$status)
)
write_tsv(enrichment_diag, "results/tables/enrichment_diagnostics.tsv")

# Main enrichment figures ------------------------------------------------------
top_n <- as.integer(cfg$reporting$enrichment_terms_each_direction)
top_terms <- function(x, n) {
  x <- x[status == "completed" & is.finite(bh_fdr) & bh_fdr < fdr_cut]
  if (!nrow(x)) return(x)
  x[order(bh_fdr, -hit_count, term_id)][seq_len(min(n, .N))]
}
f13_data <- rbind(top_terms(go_ora[ontology == "BP" & direction == "cancer_up"], top_n),
                  top_terms(go_ora[ontology == "BP" & direction == "normal_up"], top_n), fill = TRUE)
stopifnot(nrow(f13_data) > 0L)
f13_data[, display_term := paste0(term_name, " (", term_id, ")")]
f13_data[, display_term := factor(display_term, levels = rev(unique(display_term[order(direction, bh_fdr)])))]
f13 <- ggplot(f13_data, aes(gene_ratio, display_term, size = hit_count,
                            fill = -log10(pmax(bh_fdr, 1e-300)))) +
  geom_point(shape = 21, alpha = 0.95, stroke = 0.5, color = portfolio_colors[["neutral_dark"]]) +
  facet_wrap(~ direction, scales = "free_y",
             labeller = as_labeller(c(cancer_up = "Cancer-up", normal_up = "Normal-up"))) +
  scale_fill_gradient(low = portfolio_colors[["stable_light"]], high = portfolio_colors[["stable_dark"]]) +
  labs(title = "GO Biological Process over-representation", subtitle = "TREAT-priority genes; cancer-up and normal-up tested separately",
       x = "Gene ratio (hits / mapped priority genes)", y = NULL, size = "Hits", fill = expression(-log[10](BH~FDR))) +
  theme_portfolio() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_figure(f13, "F13_GO_BP_ORA_dotplot", 14, 9)
write_tsv(f13_data, "results/tables/figure_F13_source.tsv")
write_caption("F13_GO_BP_ORA_dotplot", "The universe is P1-tested Entrez genes annotated to GO BP. Point size is the hit count; ratios use the ontology-mapped priority denominator. Similar GO terms are not independent; the full unpruned table is delivered.")

gsea_pos <- gsea_dt[direction == "cancer_up" & bh_fdr < fdr_cut][order(bh_fdr, -abs(NES), term_id)][seq_len(min(top_n, .N))]
gsea_neg <- gsea_dt[direction == "normal_up" & bh_fdr < fdr_cut][order(bh_fdr, -abs(NES), term_id)][seq_len(min(top_n, .N))]
f14_data <- rbind(gsea_pos, gsea_neg, fill = TRUE)
stopifnot(nrow(f14_data) > 0L)
f14_data[, display_term := paste0(term_name, " (", term_id, ")")]
f14_data[, display_term := factor(display_term, levels = rev(unique(display_term[order(NES)])))]
f14 <- ggplot(f14_data, aes(NES, display_term, color = direction, size = -log10(pmax(bh_fdr, 1e-300)))) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.55,
             color = portfolio_colors[["neutral_mid"]]) +
  geom_point(aes(shape = direction), alpha = 0.92) +
  scale_color_manual(values = c(cancer_up = portfolio_colors[["cancer"]],
                                normal_up = portfolio_colors[["normal"]]),
                     labels = c(cancer_up = "Cancer-up", normal_up = "Normal-up")) +
  scale_shape_manual(values = c(cancer_up = 16, normal_up = 17),
                     labels = c(cancer_up = "Cancer-up", normal_up = "Normal-up")) +
  labs(title = "Reactome preranked GSEA", subtitle = "Commercially reusable fallback collection; rank is P1 moderated t",
       x = "Normalized enrichment score (NES)", y = NULL, color = "Direction", shape = "Direction",
       size = expression(-log[10](BH~FDR))) +
  theme_portfolio() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_figure(f14, "F14_Reactome_GSEA_NES", 13, 9)
write_tsv(f14_data, "results/tables/figure_F14_source.tsv")
write_caption("F14_Reactome_GSEA_NES", "Positive NES indicates enrichment toward cancer-up ranks; negative NES indicates enrichment toward normal-up ranks. Reactome is used because the MSigDB Hallmark commercial-use review was not confirmed.")

curve_pos <- gsea_dt[direction == "cancer_up" & bh_fdr < fdr_cut][order(bh_fdr, -abs(NES), term_id)][seq_len(min(3L, .N))]
curve_neg <- gsea_dt[direction == "normal_up" & bh_fdr < fdr_cut][order(bh_fdr, -abs(NES), term_id)][seq_len(min(3L, .N))]
curve_terms <- rbind(curve_pos, curve_neg, fill = TRUE)
stopifnot(nrow(curve_terms) > 0L)
curve_plots <- lapply(seq_len(nrow(curve_terms)), function(i) {
  id <- curve_terms$term_id[i]
  direction_color <- if (curve_terms$direction[i] == "cancer_up") {
    portfolio_colors[["cancer_dark"]]
  } else {
    portfolio_colors[["normal_dark"]]
  }
  p <- plotEnrichment(rx_pathways[[id]], ranks, gseaParam = gsea_exp)
  p$layers[[1]]$aes_params$colour <- direction_color
  p$layers[[1]]$aes_params$linewidth <- 0.8
  p$layers[[2]]$aes_params$colour <- portfolio_colors[["neutral_mid"]]
  p$layers[[2]]$aes_params$alpha <- 0.55
  p$layers[[3]]$aes_params$colour <- direction_color
  p$layers[[4]]$aes_params$colour <- direction_color
  p$layers[[3]]$aes_params$linewidth <- 0.45
  p$layers[[4]]$aes_params$linewidth <- 0.45
  p$layers[[3]]$aes_params$linetype <- "22"
  p$layers[[4]]$aes_params$linetype <- "22"
  p$layers[[5]]$aes_params$colour <- portfolio_colors[["neutral_dark"]]
  p$layers[[5]]$aes_params$linewidth <- 0.45
  p +
    labs(title = curve_terms$term_name[i], subtitle = sprintf("%s | NES %.2f | BH FDR %.2g", id, curve_terms$NES[i], curve_terms$bh_fdr[i]),
         x = "Rank in P1 moderated-t vector", y = "Running enrichment score") + theme_portfolio() +
    theme(legend.position = "none", plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 8),
          plot.margin = margin(7, 9, 7, 7))
})
curve_grob <- arrangeGrob(grobs = curve_plots, ncol = 2)
draw_f15 <- function() { grid::grid.newpage(); grid::grid.draw(curve_grob) }
save_base_figure("F15_selected_enrichment_curves", draw_f15, 13, 10)
curve_source <- rbindlist(lapply(seq_len(nrow(curve_terms)), function(i) {
  id <- curve_terms$term_id[i]
  data.table(term_id = id, term_name = curve_terms$term_name[i], rank = seq_along(ranks),
             entrez_id = names(ranks), moderated_t = unname(ranks), in_gene_set = names(ranks) %in% rx_pathways[[id]])
}))
write_tsv(curve_source, "results/tables/figure_F15_source.tsv")
write_caption("F15_selected_enrichment_curves", "Up to three BH-significant Reactome sets in each direction are selected by FDR, absolute NES, then stable pathway ID. Curves use all ranked P1-tested genes; leading edges are listed in the full GSEA table.")

message(sprintf("[04] PASS: GO BP %d significant rows; Reactome ORA %d; Reactome GSEA %d; rank ties %.4f%%",
                go_ora[ontology == "BP" & bh_fdr < fdr_cut, .N], rx_ora[bh_fdr < fdr_cut, .N],
                gsea_dt[bh_fdr < fdr_cut, .N], 100 * tie_fraction))

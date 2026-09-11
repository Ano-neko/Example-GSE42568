# Clean rebuild verification

**Status: PASS (10/10 key outputs byte-identical)**  
**Rebuild completed:** 2026-09-11 08:29:19 UTC  
**Runner:** `bash scripts/run_all.sh` from the repository root

Generated metadata, derived objects, results, reports and output directories were moved to a recoverable temporary backup. The pipeline then rebuilt from the retained checksum-locked official GEO matrix and GPL570 annotation input. All stages, 13/13 reproducibility guards and 32/32 Scientific QA checks passed.

The following machine-readable key outputs were SHA-256 compared to their pre-rebuild versions:

- `data/metadata/sample_manifest.tsv` — identical
- `results/tables/matrix_integrity.tsv` — identical
- `results/tables/probe_selection_audit.tsv` — identical
- `results/tables/gene_annotation.tsv` — identical
- `results/tables/de_primary_all_genes.tsv` — identical
- `results/tables/de_priority_genes.tsv` — identical
- `results/tables/de_sensitivity_summary.tsv` — identical
- `results/tables/go_ora_all_terms.tsv` — identical (retained locally; the content-equivalent `.tsv.gz` is the GitHub copy)
- `results/tables/pathway_ora_all_terms.tsv` — identical
- `results/tables/gsea_all_gene_sets.tsv` — identical

The complete hashes are in `results/tables/rebuild_verification.tsv`. Report PDFs were regenerated and validated separately for page count, extractable text and rendered-page layout; PDF bytes are not used as deterministic statistical-result checks.

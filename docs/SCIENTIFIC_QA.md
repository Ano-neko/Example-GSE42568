# Scientific QA

**Overall status: PASS_WITH_DOCUMENTED_LIMITATIONS**  
**Review date:** 2026-09-11  
**Contrast:** cancer minus normal breast tissue

## Outcome

- 32 checks evaluated; 32 passed; 0 hard failures; 0 soft failures.
- P1 tested 20357 genes; 5381 moderated BH discoveries; 2835 TREAT priorities; 2745 headline-eligible genes.
- Design stability: Green. P1/S1 Spearman 0.936; P1/S2 Spearman 1.000; both preserve 100% of P1-priority directions.
- S3 is formally not applicable because no sample met the signed technical-failure definition. Four reviewed QC candidates were retained in P1 and jointly removed only in a separately labelled influence sensitivity.
- KEGG and MSigDB Hallmark were not run because the required commercial-use licence review was not documented. Reactome is the licensed pathway/GSEA fallback.

## Interpretation boundary

This analysis estimates an adjusted cross-sectional association in heterogeneous bulk tissues. It does not establish causality, diagnosis, prognosis, treatment response, pathway activation, or cell-intrinsic regulation. Normal-tissue provenance and age are unavailable; run is partially confounded with group; and processed GEO data cannot support raw-array NUSE/RLE, degradation, or image QC.

## Machine-readable evidence

- `results/tables/acceptance_checks.tsv`
- `results/tables/biological_sanity_checks.tsv`
- `results/tables/de_sensitivity_summary.tsv`
- `results/tables/leave_one_run_out_summary.tsv`
- `results/tables/gene_set_manifest.tsv`


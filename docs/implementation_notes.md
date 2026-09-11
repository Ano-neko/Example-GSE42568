# Implementation notes and bounded omissions

## Completed core chain

Official input validation, sample manifest, processed-data QC, signed sample review, current probe annotation, representative and median gene matrices, P1/S1/S2/S4 models, QC-candidate influence analysis, probe-level supplement, leave-one-run-out analysis, GO ORA, Reactome ORA/GSEA, F01–F15, Scientific QA, HTML/PDF reports, environment locks, artifact manifest and clean-rebuild verification are complete.

The specification-named `results/tables/go_ora_all_terms.tsv` is retained in the local delivery. It is larger than GitHub's 100 MB per-file limit, so the pipeline also writes a content-equivalent `go_ora_all_terms.tsv.gz`; only the compressed copy is intended for version control.

Reproducibility is provided by the pinned Conda environment in `config/environment.yml` and `config/conda-explicit-linux-64.txt`. Runtime commands default to `Rscript`, `python` and `conda` on the active `PATH`; `BIOINFO_ENV_PREFIX` is an optional explicit override. No project path or environment prefix is fixed in executable code, and no separate renv lock is maintained.

## Deliberately bounded or omitted

- Family SOFT was not retained because all fields required by v1 are available in the checksum-locked series-matrix header. This is recorded in `source_manifest.tsv`; it did not substitute another biological input.
- Formal S3 is not applicable because no sample met the signed technical-failure definition. Four reviewed candidates are handled in a separately named influence sensitivity.
- KEGG and MSigDB Hallmark were not accessed because commercial-use permission/review was not documented. Status rows are delivered; Reactome is the named fallback.
- Full-data CI is not run on a hosted service. Local guards test the synthetic contrast, schemas, annotation exclusions, universe membership and licence gates; the controlled runner performs the full rebuild.
- `scripts/run_all.sh` is the acceptance runner. `_targets.R` is supplied as a DAG interface, but the verified v1 release used the transparent sequential runner to avoid making workflow-engine refinement a release blocker.

## Implementation issue resolved before release

An early Scientific QA attempt used two incorrect checker expressions: stale matrix metric keys and a prohibited-module search that matched its own rule text. Neither affected inputs, fitted models or result tables. Both checks were corrected; the final run passes 32/32 and the clean rebuild reproduces all 10 key result tables byte-for-byte.

# QC candidate review

**Status:** reviewed and released for primary analysis  
**Reviewer:** Scientific Architect  
**Review date:** 2026-09-11

## Trigger

Four samples triggered at least two independent processed-expression warning categories. This was a pre-specified pause/review condition, not an automatic exclusion rule.

## Evidence considered

- All four arrays have smooth, finite expression distributions within the continuous range of the study; none shows a grossly truncated, singular, or globally corrupted profile.
- Pearson connectivity and robust PCA confirm that the samples are unusual, especially GSM1045193, GSM1045217, and GSM1045302, but each retains plausible biological nearest neighbours.
- GSM1045192 is supported by distribution and PCA warnings but does not cross the correlation warning boundary.
- The deposited matrix provides no raw CEL quality metrics (for example RNA degradation, NUSE/RLE, or image artefacts) that could independently establish technical failure.
- The warning pattern is not sufficient to distinguish technical quality from true tissue heterogeneity, particularly because cancer and normal breast tissue differ substantially and run is partially confounded with tissue group.

## Decision

All 121 samples are retained in the primary model. GSM1045192, GSM1045193, GSM1045217, and GSM1045302 are excluded jointly in the locked QC-candidate sensitivity model. This deliberately tests influence without retroactively labelling any sample as a failed array.

A primary conclusion is considered robust only if its direction and magnitude remain materially consistent in this and the other pre-specified sensitivity analyses. The machine-readable signed decision is in `results/tables/sample_exclusion_ledger.tsv`; the supporting values are in `results/tables/qc_candidate_review.tsv`.


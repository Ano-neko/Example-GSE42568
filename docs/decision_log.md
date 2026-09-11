# Decision log

## 2026-09-11 — v1 locked implementation

- Accepted the checksum-locked GEO processed matrix; no second log or normalization.
- Locked P1 to run proxy plus tissue group, normal reference, cancer-minus-normal coefficient.
- Used current Entrez annotation and an outcome-blind representative-probe rule.
- Retained four reviewed QC candidates in P1; no raw evidence established failure. Removed them only in an influence model; formal S3 is not applicable.
- Completed group-only, mixed-run, no-trend, group-aware, median-collapse, probe-level and leave-one-run-out sensitivities; final traffic light Green.
- Kept TREAT 1.2-fold priority and 1.1/1.5 count sensitivities unchanged.
- Ran GO and Reactome with resource-specific tested-and-annotated universes.
- Did not access KEGG or MSigDB content because commercial-use review was not documented.
- Kept family SOFT as a documented infrastructure omission; scope was not expanded.

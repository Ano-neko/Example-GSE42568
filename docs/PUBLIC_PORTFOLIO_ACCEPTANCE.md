# Public portfolio edition acceptance

## Release logic

**Full analysis run passed 32/32 QA.**

**Public portfolio edition intentionally omits large/generated artifacts.**

**The included public edition passed its own package-integrity check: 7/7 checks passed.**

The full-workspace Scientific QA and the public-edition acceptance answer different questions. QA031 proves that figure source tables existed and were audited in the complete analysis workspace. The public edition then intentionally excludes those regenerable tables, all `data/`, run-specific logs, the uncompressed GO duplicate, and full-workspace manifest/checksum files.

The retained public edition contains 147 files under the declared inclusion policy. It preserves code, frozen Conda specifications, tests, documentation, client reports, PNG/PDF figures with captions, compact/key statistical tables, licence records and both full-workspace and public-edition QA ledgers.

## Package-integrity checks

| Check | Requirement | Status | Evidence |
|---|---|---|---|
| PKG001 | Full analysis QA passed before public-edition filtering | PASS | 32/32 full-workspace checks passed |
| PKG002 | QA031 is explicitly scoped to the complete analysis workspace | PASS | Complete analysis workspace contains all F01-F15 PDF/PNG/caption and source-data files |
| PKG003 | Required public code, reports and compact results exist | PASS | No required files missing |
| PKG004 | Intentional public-edition omissions are declared | PASS | All omission rules present |
| PKG005 | Public-edition text files contain no server-specific project or environment path | PASS | No machine-specific path found |
| PKG006 | Pinned Conda records are the sole environment lock | PASS | renv.lock absent; Conda source of truth documented |
| PKG007 | Public edition includes PNG, PDF and caption for all 16 released figures | PASS | 48/48 expected public figure files present and nonempty |

Machine-readable evidence: `results/tables/public_portfolio_acceptance.tsv`.

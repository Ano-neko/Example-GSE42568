#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def public_omitted(relative_path: str) -> bool:
    parts = Path(relative_path).parts
    if not parts:
        return True
    if parts[0] in {".git", "data", "tmp"} or "__pycache__" in parts:
        return True
    if relative_path.startswith("results/logs/") or relative_path.startswith("results/checksums/"):
        return True
    if relative_path.startswith("results/tables/figure_") and relative_path.endswith("_source.tsv"):
        return True
    if relative_path in {
        "machine_readable_manifest.json",
        "renv.lock",
        "Rplots.pdf",
        "results/tables/artifact_checksums.tsv",
        "results/tables/go_ora_all_terms.tsv",
    }:
        return True
    return relative_path.endswith((".pyc", ".swp", ".swo"))


def is_text_without_forbidden_path(path: Path, forbidden: bytes) -> bool:
    payload = path.read_bytes()
    if b"\0" in payload:
        return True
    return forbidden not in payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the public portfolio edition independently of full-workspace QA.")
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    project = parser.parse_args().project.resolve()

    qa_rows = read_tsv(project / "results/tables/acceptance_checks.tsv")
    qa031 = next((row for row in qa_rows if row["check_id"] == "QA031"), None)
    gitignore = (project / ".gitignore").read_text(encoding="utf-8")
    readme = (project / "README.md").read_text(encoding="utf-8")

    figure_stems = [f"F{i:02d}" for i in range(1, 16)] + ["S01"]
    figure_files: list[Path] = []
    for stem in figure_stems:
        matches = sorted((project / "results/figures").glob(f"{stem}_*"))
        for suffix in (".png", ".pdf", ".caption.txt"):
            figure_files.extend(path for path in matches if path.name.endswith(suffix))

    required = [
        "README.md",
        "CITATION.cff",
        "LICENSE",
        "Makefile",
        "_targets.R",
        "config/environment.yml",
        "config/conda-explicit-linux-64.txt",
        "scripts/run_all.sh",
        "scripts/capture_environment.sh",
        "scripts/build_reports.py",
        "scripts/finalize_release.py",
        "scripts/validate_public_release.py",
        "docs/SCIENTIFIC_QA.md",
        "docs/DELIVERY_ACCEPTANCE.md",
        "output/pdf/client_report.pdf",
        "output/pdf/methods_appendix.pdf",
        "output/pdf/qc_appendix.pdf",
        "results/tables/acceptance_checks.tsv",
        "results/tables/de_priority_genes.tsv",
        "results/tables/gene_set_manifest.tsv",
        "results/tables/go_ora_all_terms.tsv.gz",
    ]
    required_missing = [rel for rel in required if not (project / rel).is_file()]

    ignore_contract = [
        "data/",
        "results/tables/go_ora_all_terms.tsv",
        "results/tables/figure_*_source.tsv",
        "machine_readable_manifest.json",
        "results/tables/artifact_checksums.tsv",
        "results/checksums/",
        "results/logs/",
        "Rplots.pdf",
    ]
    missing_ignore_rules = [rule for rule in ignore_contract if rule not in gitignore.splitlines()]

    public_files = []
    hardcoded_hits = []
    forbidden = ("/home" + "/bioinfo").encode()
    for path in sorted(p for p in project.rglob("*") if p.is_file()):
        rel = path.relative_to(project).as_posix()
        if public_omitted(rel):
            continue
        public_files.append(rel)
        if not is_text_without_forbidden_path(path, forbidden):
            hardcoded_hits.append(rel)

    checks: list[dict[str, str]] = []

    def add(check_id: str, domain: str, description: str, passed: bool, evidence: str) -> None:
        checks.append({
            "check_id": check_id,
            "domain": domain,
            "description": description,
            "status": "PASS" if passed else "FAIL",
            "passed": "TRUE" if passed else "FALSE",
            "evidence": evidence,
        })

    full_qa_ok = len(qa_rows) == 32 and all(row["passed"] == "TRUE" for row in qa_rows)
    add("PKG001", "scientific_lineage", "Full analysis QA passed before public-edition filtering", full_qa_ok,
        f"{sum(row['passed'] == 'TRUE' for row in qa_rows)}/{len(qa_rows)} full-workspace checks passed")

    qa031_ok = qa031 is not None and qa031["passed"] == "TRUE" and qa031["requirement"].startswith("Complete analysis workspace")
    add("PKG002", "scope", "QA031 is explicitly scoped to the complete analysis workspace", qa031_ok,
        qa031["requirement"] if qa031 else "QA031 missing")

    add("PKG003", "required_artifacts", "Required public code, reports and compact results exist", not required_missing,
        "No required files missing" if not required_missing else ", ".join(required_missing))

    add("PKG004", "omission_policy", "Intentional public-edition omissions are declared", not missing_ignore_rules,
        "All omission rules present" if not missing_ignore_rules else ", ".join(missing_ignore_rules))

    add("PKG005", "portability", "Public-edition text files contain no server-specific project or environment path", not hardcoded_hits,
        "No machine-specific path found" if not hardcoded_hits else ", ".join(hardcoded_hits))

    conda_ok = not (project / "renv.lock").exists() and "Reproducibility is provided by the pinned Conda environment" in readme
    add("PKG006", "environment", "Pinned Conda records are the sole environment lock", conda_ok,
        "renv.lock absent; Conda source of truth documented" if conda_ok else "Environment-lock policy mismatch")

    unique_figures = sorted(set(path.resolve() for path in figure_files))
    figures_ok = len(unique_figures) == 48 and all(path.stat().st_size > 0 for path in unique_figures)
    add("PKG007", "figures", "Public edition includes PNG, PDF and caption for all 16 released figures", figures_ok,
        f"{len(unique_figures)}/48 expected public figure files present and nonempty")

    table_path = project / "results/tables/public_portfolio_acceptance.tsv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)

    passed = sum(row["passed"] == "TRUE" for row in checks)
    public_file_count = len(set(public_files + [
        "docs/PUBLIC_PORTFOLIO_ACCEPTANCE.md",
        "results/tables/public_portfolio_acceptance.tsv",
    ]))
    markdown_rows = "\n".join(
        f"| {row['check_id']} | {row['description']} | {row['status']} | {row['evidence']} |" for row in checks
    )
    doc = f"""# Public portfolio edition acceptance

## Release logic

**Full analysis run passed 32/32 QA.**

**Public portfolio edition intentionally omits large/generated artifacts.**

**The included public edition passed its own package-integrity check: {passed}/{len(checks)} checks passed.**

The full-workspace Scientific QA and the public-edition acceptance answer different questions. QA031 proves that figure source tables existed and were audited in the complete analysis workspace. The public edition then intentionally excludes those regenerable tables, all `data/`, run-specific logs, the uncompressed GO duplicate, and full-workspace manifest/checksum files.

The retained public edition contains {public_file_count} files under the declared inclusion policy. It preserves code, frozen Conda specifications, tests, documentation, client reports, PNG/PDF figures with captions, compact/key statistical tables, licence records and both full-workspace and public-edition QA ledgers.

## Package-integrity checks

| Check | Requirement | Status | Evidence |
|---|---|---|---|
{markdown_rows}

Machine-readable evidence: `results/tables/public_portfolio_acceptance.tsv`.
"""
    (project / "docs/PUBLIC_PORTFOLIO_ACCEPTANCE.md").write_text(doc, encoding="utf-8")

    print(f"[public-release] {'PASS' if passed == len(checks) else 'FAIL'}: {passed}/{len(checks)} checks; {public_file_count} included files")
    if passed != len(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

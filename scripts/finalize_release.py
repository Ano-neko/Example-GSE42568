#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import mimetypes
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def category(rel: str) -> str:
    if rel.startswith("results/figures/"):
        return "figure"
    if rel.startswith("results/tables/"):
        return "table"
    if rel.startswith("output/pdf/") or rel.startswith("reports/"):
        return "report"
    if rel.startswith("R/") or rel.startswith("scripts/") or rel.startswith("tests/"):
        return "executable_source"
    if rel.startswith("config/"):
        return "configuration_environment"
    if rel.startswith("docs/"):
        return "documentation"
    return "project_metadata"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    project = parser.parse_args().project.resolve()
    excluded_roots = {".git", "tmp"}
    excluded_prefixes = ("data/raw/", "data/derived/", "results/logs/")
    excluded_names = {"machine_readable_manifest.json", "artifact_checksums.tsv", "artifact_checksums.sha256"}

    def included(path: Path) -> bool:
        rel = path.relative_to(project).as_posix()
        if any(part in excluded_roots or part == "__pycache__" for part in path.relative_to(project).parts):
            return False
        if rel.startswith(excluded_prefixes) or path.name in excluded_names or path.suffix == ".pyc":
            return False
        return path.is_file()

    de = read_tsv(project / "results/tables/de_model_diagnostics.tsv")
    de_map = {r["metric"]: r["value"] for r in de}
    enrichment = read_tsv(project / "results/tables/enrichment_diagnostics.tsv")
    enrichment_map = {r["metric"]: r["value"] for r in enrichment}
    artifacts = []
    for path in sorted((p for p in project.rglob("*") if included(p)), key=lambda p: p.relative_to(project).as_posix()):
        rel = path.relative_to(project).as_posix()
        artifacts.append({
            "path": rel,
            "size_bytes": path.stat().st_size,
            "sha256": digest(path),
            "media_type": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
            "category": category(rel),
        })

    manifest = {
        "schema_version": "1.0.0",
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "project": {"accession": "GSE42568", "platform": "GPL570", "release": "v1.0",
                    "contrast": "cancer-minus-normal", "scientific_status": "PASS_WITH_DOCUMENTED_LIMITATIONS"},
        "input": {"matrix_sha256": "43b93f72f1d81838dcd2a5b6ccc3bc1875e46230dad0d7346236773a99e32b76",
                  "samples": 121, "cancer": 104, "normal": 17, "probe_sets": 54675},
        "results": {"tested_genes": int(float(de_map["tested_genes"])),
                    "moderated_bh_discoveries": int(float(de_map["statistically_discovered_genes"])),
                    "treat_1_2_priority_genes": int(float(de_map["priority_genes"])),
                    "headline_eligible_genes": int(float(de_map["headline_eligible_genes"])),
                    "design_traffic_light": de_map["final_design_traffic_light"],
                    "go_bp_significant_rows": int(float(enrichment_map["go_bp_significant_rows"])),
                    "reactome_gsea_significant_rows": int(float(enrichment_map["reactome_gsea_significant_rows"]))},
        "licensing": {"GO": "completed", "Reactome": "completed",
                      "KEGG": enrichment_map["kegg_status"], "MSigDB_Hallmark": enrichment_map["msigdb_hallmark_status"]},
        "checksum_scope_note": "Raw/derived data, run-specific logs, tmp, git internals, and checksum files are excluded. This is the full-workspace manifest; the public portfolio edition has a separate acceptance ledger.",
        "artifacts": artifacts,
    }
    manifest_path = project / "machine_readable_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # Build the final checksum ledger including this manifest but excluding the checksum files themselves.
    checksum_paths = sorted(
        [p for p in project.rglob("*") if included(p)] + [manifest_path],
        key=lambda p: p.relative_to(project).as_posix(),
    )
    rows = [(p.relative_to(project).as_posix(), p.stat().st_size, digest(p)) for p in checksum_paths]
    table_path = project / "results/tables/artifact_checksums.tsv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        writer.writerows(rows)
    sha_path = project / "results/checksums/artifact_checksums.sha256"
    sha_path.parent.mkdir(parents=True, exist_ok=True)
    sha_path.write_text("".join(f"{sha}  {rel}\n" for rel, _, sha in rows), encoding="utf-8")
    print(f"[finalize] manifest artifacts={len(artifacts)}; checksum rows={len(rows)}")


if __name__ == "__main__":
    main()

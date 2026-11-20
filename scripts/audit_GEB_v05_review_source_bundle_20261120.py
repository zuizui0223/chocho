#!/usr/bin/env python3
"""Fail-closed source-lineage and reviewer-archive feasibility audit for GEB v0.5.

This script *does not* produce a peer-review-ready or author-anonymous archive.
It checks the actual four-figure paper against its pinned source commit and
builds a deliberately labelled INTERNAL derived-evidence preview, preserving
source identity for scientific reproducibility. It never claims raw occurrence
or historical transport data exist just because aggregate receipts exist.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

SOURCE_COMMIT = "3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"
MANUSCRIPT = Path("manuscript/integrated_butterfly_resource_access_20261109")
FIGURE_DIR = Path("figures/integrated_island_port_butterfly_v03")
FIGURE_NUMBERS = FIGURE_DIR / "FIGURE_V04_SOURCE_NUMBERS.json"
MAIN_FIGURES = (
    "Fig1_resource_geography_scales_v03.svg",
    "Fig2_ports_native_alien_floras_isolation_v05.svg",
    "Fig3_global_host_contributors_at_ports_restructured_v04.svg",
    "Fig4_independent_butterfly_record_predictability_v04.svg",
)
SI_FIGURES = (
    "Fig2_global_contributors_at_ports.svg",
    "Fig4_independent_butterfly_record_predictability_full_south_america_v04.svg",
)
PAPER_PATHS = [
    MANUSCRIPT / "ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md",
    MANUSCRIPT / "GEB_SUPPORTING_INFORMATION_V05.md",
    MANUSCRIPT / "FIGURE_CAPTIONS_V03.md",
    MANUSCRIPT / "GEB_V05_DATA_AND_SUBMISSION_HANDOFF_20261120.md",
    FIGURE_NUMBERS,
    FIGURE_DIR / "FIGURE4_MARINE_SOURCE_RECEIPT.json",
    Path("results/island-mainland-resource-q1/region_metrics.csv"),
    Path("results/island-mainland-resource-q1/island_mainland_summary.json"),
    Path("results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv"),
    Path("results/isolation_port_resource_reanalysis_20261120/SOURCE239_PAPER_RECEIPT_V01.json"),
    Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv"),
]
PAPER_PATHS.extend(FIGURE_DIR / p for p in MAIN_FIGURES + SI_FIGURES)
CODE_PATHS = [
    Path("scripts/analyze_butterfly_anthropogenic_resource_expansion.py"),
    Path("scripts/analyze_butterfly_host_contribution_concentration.py"),
    Path("scripts/analyze_butterfly_host_breadth_geography.py"),
    Path("scripts/analyze_gbif_global224_observed_beta.py"),
    Path("scripts/aggregate_gbif_global224_wgsrpd3_annual.py"),
    Path("scripts/audit_island_port_isolation_239_20261120.py"),
    Path("scripts/audit_island_port_spatial_239_20261120.py"),
    Path("scripts/verify_integrated_geb_v04_editorial.py"),
    Path("scripts/build_GEB_v05_blinded_main_and_SI.py"),
]
# Supplemental source evidence intentionally kept in the pinned-source checkout,
# rather than mistakenly searched for in the newer journal paper branch.
FROZEN_SOURCE_PATHS = [
    Path("provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json"),
]
NOT_IN_RELEASE = [
    "83,573-row original butterfly-island-year Parquet and source record IDs",
    "225-butterfly × 51-landmass original GBIF occurrence datasetKey attribution",
    "1,596-port botanical occurrence-level GBIF datasetKey attribution",
    "original WPI 2019 input file and independently verified original-to-GeoJSON conversion",
    "full original records for independently dated 155-species retrospective prediction",
    "reviewer-visible anonymous archive location / peer-review private link",
    "corresponding author-approved identifying title page, disclosures and media permissions",
]
IDENTITY_PATTERNS = {
    "repository_account": re.compile(rb"zuizui0223", re.I),
    "email_address": re.compile(rb"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}", re.I),
    "home_directory": re.compile(rb"/home/[^ /]+/", re.I),
    "identifying_credit": re.compile(rb"ZHANG[ _-]?(?:RUIQI|RACHEL)", re.I),
}
# No source account name or correspondence identity should be silently hidden
# by rewrites that would invalidate source provenance. Scan and block release.
def digest(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def status_of(root: Path, path: Path, group: str) -> dict:
    p = root / path
    if not p.is_file():
        return {"origin": group, "path": str(path), "exists": False}
    data = p.read_bytes()
    return {"origin": group, "path": str(path), "exists": True,
            "size_bytes": len(data), "sha256": digest(data)}


def read_json(root: Path, path: Path) -> dict:
    return json.loads((root / path).read_text(encoding="utf-8"))


def canonical_zip(name: str, data: bytes, archive: ZipFile) -> None:
    info = ZipInfo(name, (2026, 11, 20, 0, 0, 0))
    info.compress_type = ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    archive.writestr(info, data)


def audit(current: Path, source: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    article = (current / PAPER_PATHS[0]).read_text(encoding="utf-8")
    si = (current / PAPER_PATHS[1]).read_text(encoding="utf-8")
    nums = read_json(current, FIGURE_NUMBERS)
    original = nums.get("main_evidence_commit_at_read")
    if original != SOURCE_COMMIT:
        raise RuntimeError(f"Figure source pin changed unexpectedly: {original}")
    if len(re.findall(r"^## Appendix S[1-8]\.", si, re.M)) != 8:
        raise RuntimeError("Current GEB S1-S8 supporting information structure drifted")
    if "### S3.3. Additional spatial autocorrelation" not in si:
        raise RuntimeError("Current GEB S3.3 Moran/MEM10 source check missing")
    if "### 3.2." not in article or "Figure 2D" not in article:
        raise RuntimeError("Current paper missing isolation panel D")
    if not all(m in article for m in ("+0.425", "+0.395", "p=0.068")):
        raise RuntimeError("Current source missing paired geographic caveats")
    counts = nums["figures"]["fig2"]["isolation_panel"]["cohort"]
    if (counts["islands"], counts["fixed_butterflies"], counts["active_islands"]) != (126, 239, 123):
        raise RuntimeError("Fixed-239 source cohort drift")
    b = nums["figures"]["fig1"]
    if (b["relative_increase_percent"]["island"], b["relative_increase_percent"]["mainland"]) != (118.9, 40.9):
        raise RuntimeError("Original island resource baseline drift")
    original_paths = [Path(s) for s in nums["source_paths"]]
    inventory = []
    for p in PAPER_PATHS + CODE_PATHS:
        inventory.append(status_of(current, p, "current_paper"))
    for p in original_paths + FROZEN_SOURCE_PATHS:
        inventory.append(status_of(source, p, "immutable_numeric_source"))
    absent = [q for q in inventory if not q["exists"]]
    if absent:
        raise RuntimeError("Missing directly required paper/source receipts: " +
                           json.dumps(absent, ensure_ascii=False))
    # The original source receipt is intentionally from the immutable checkout,
    # not a similarly named but unrelated post-hoc paper-branch file.
    for name in MAIN_FIGURES + SI_FIGURES:
        ET.fromstring((current / FIGURE_DIR / name).read_bytes())
    staged = [(("current_paper", p), (current / p).read_bytes())
              for p in PAPER_PATHS + CODE_PATHS]
    staged += [(("immutable_numeric_source", p), (source / p).read_bytes())
               for p in original_paths + FROZEN_SOURCE_PATHS]
    leaks = []
    for (origin, p), data in staged:
        if p.suffix.lower() not in (".py", ".md", ".json", ".csv", ".svg"):
            continue
        names = [kind for kind, pattern in IDENTITY_PATTERNS.items()
                 if pattern.search(data)]
        if names:
            leaks.append({"origin": origin, "file": str(p), "identity_categories": names})
    # Additional explicit block: 8 immutable figure summaries and selected
    # aggregate data are not the original record-level holdings for all panels.
    gates = {
        "main_and_si_visual_evidence": "PASS_FOR_INTERNAL_QA",
        "source_figures_reconciled_to_immutable_snapshot": "PASS_FOR_INTERNAL_QA",
        "candidate_text_and_scripts_checked_for_identity": "SCAN_COMPLETE_NOT_ANONYMIZED",
        "all_original_data_and_record_ids": "BLOCKED",
        "dataset_publisher_attribution_for_all_GBIF_queries": "BLOCKED",
        "WPI_original_conversion_and_rights": "BLOCKED",
        "anonymous_reviewer_access_and_stable_archive": "BLOCKED",
        "coauthor_and_corresponding_author_approval": "BLOCKED"
    }
    report = {
        "schema": "chocho_GEB_v05_review_source_preflight_20261120_v1",
        "state": "INTERNAL_EVIDENCE_PREVIEW_ONLY_NOT_JOURNAL_REVIEW_ARCHIVE",
        "paper_ref_as_checked_out": "current",
        "pinned_numeric_source_git_commit": SOURCE_COMMIT,
        "scientific_cohort": counts,
        "figures": {"main_svg_count": len(MAIN_FIGURES), "SI_svg_count": len(SI_FIGURES),
                    "source_numeric_receipts": len(original_paths)},
        "included_inventory_count": len(inventory),
        "included_inventory": inventory,
        "identity_scanner": {"candidate_files_with_identifiers": len(leaks),
                             "hits": leaks},
        "unresolved_inputs_or_approvals": NOT_IN_RELEASE,
        "publication_gates": gates,
        "not_represented": [
            "third-party raw WCVP, HOSTS and GBIF occurrence records",
            "historical WPI original conversion, vessel arrivals or individual plant introduction dates",
            "butterfly larval survival, population gain or observed island establishment",
        ],
    }
    with (out / "SOURCE_PREVIEW_MANIFEST.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=("origin","path","exists","size_bytes","sha256"))
        writer.writeheader()
        for row in inventory:
            writer.writerow(row)
    (out / "GEB_V05_SOURCE_PREVIEW_AUDIT.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    guide = """# GEB v0.5 internal source evidence preview — NOT an anonymous peer-review deposit

This is a limited, **non-anonymized** source-evidence preview generated by an
automated source-pinned audit. It contains selected code, figure SVGs, numeric
receipts and aggregate island resource tables; it intentionally omits original
record-level GBIF, HOSTS/WCVP, private historical Parquet, unreconciled WPI
source conversion and nonlicensed third-party media.

Not every analysis in the GEB paper can be rerun from this ZIP alone. Do not
upload this ZIP to the Wiley editorial system or attach it to an anonymous
reviewer data statement. For missing source families, see the sibling JSON
audit's 'unresolved_inputs_or_approvals'. The exact published source-version
attributions are HOSTS DOI 10.5519/havt50xw, WCVP v13 DOI
10.34885/0yex-xv26 and separate observed 224-species GBIF DOI
10.15468/dl.umqmsf; the last is **not** a DOI for 225×51 or 1,596-port
GBIF acquisitions. Current research article four main figures are separate
from two SI figures and from the source-manuscript v0.2 submission scripts.

Source attribution and data licensing still require human confirmation before
any sharing beyond this internal GitHub Actions QA. No DOI was minted.
"""
    reportname = "README_INTERNAL_NOT_FOR_ANONYMOUS_REVIEW.md"
    with ZipFile(out / "GEB_V05_INTERNAL_SOURCE_PREVIEW_NOT_ANONYMOUS.zip", "w") as z:
        canonical_zip(reportname, guide.encode("utf-8"), z)
        for (origin, p), data in sorted(staged, key=lambda item: (item[0][0], str(item[0][1]))):
            canonical_zip(f"{origin}/{p.as_posix()}", data, z)
        canonical_zip("GEB_V05_SOURCE_PREVIEW_AUDIT.json",
                      (out / "GEB_V05_SOURCE_PREVIEW_AUDIT.json").read_bytes(), z)
        canonical_zip("SOURCE_PREVIEW_MANIFEST.csv",
                      (out / "SOURCE_PREVIEW_MANIFEST.csv").read_bytes(), z)
    with ZipFile(out / "GEB_V05_INTERNAL_SOURCE_PREVIEW_NOT_ANONYMOUS.zip") as z:
        assert len(z.namelist()) == len(staged) + 3
        assert all(name.startswith(("current_paper/", "immutable_numeric_source/"))
                   or name in (reportname, "GEB_V05_SOURCE_PREVIEW_AUDIT.json",
                               "SOURCE_PREVIEW_MANIFEST.csv") for name in z.namelist())
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--current-root", required=True, type=Path)
    ap.add_argument("--source-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    report = audit(args.current_root, args.source_root, args.out_dir)
    print("GEB_V05_INTERNAL_SOURCE_PREVIEW_QA_PASS")
    for hit in report["identity_scanner"]["hits"]:
        print("REVIEW_IDENTITY_LEAK_FILE", hit["origin"], hit["file"], ",".join(hit["identity_categories"]))
    print(json.dumps({
        "state": report["state"],
        "inventory_count": report["included_inventory_count"],
        "main_and_SI_figures": report["figures"],
        "identity_leak_files": report["identity_scanner"]["candidate_files_with_identifiers"],
        "blocked_gates": [k for k, v in report["publication_gates"].items()
                          if v == "BLOCKED"]}, ensure_ascii=False))

if __name__ == "__main__":
    main()

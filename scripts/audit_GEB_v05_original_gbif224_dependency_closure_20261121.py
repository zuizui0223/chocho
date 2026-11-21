#!/usr/bin/env python3
"""Source-locked AST audit of the ORIGINAL 224-species GBIF model code lineage.

No original occurrence records are downloaded or published here. All source
modules are located through exact archival checksums, not reimplemented.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path

ROOT_SCRIPTS = (
    "aggregate_gbif_global224_wgsrpd3_annual.py",
    "analyze_gbif_global224_observed_beta.py",
)
# Four directly imported missing modules and their two transitive imports were
# located in the unmerged historical source branch, fixed by immutable commit.
SOURCE_COMMIT = "dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0"
PAPER_COMMIT_ORIGINAL = "3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"
REQUIRED_SOURCE_BLOBS = {
    "analyze_documented_host_specificity.py": "f024f1a85ff95b2a0ba2c43fb5e2cac6ee075510",
    "analyze_matched_glonaf_consumer_history.py": "95f274ecda2e7ae1e85259d0070ad9920fe9a61b",
    "analyze_colonial_legacy_resource_beta.py": "1bfdb98ed331c64cc72bcc84f2200b176d229b9e",
    "analyze_gbif_31_observed_host_beta.py": "f6dbb8407ca5c3bfdaf64339bf92a1ce68a9d71f",
    "build_lenzner_tdwg_colonial_crosswalk.py": "f9ac05d73f4a1de8ca27fa14fe2e6d790b411c34",
    "audit_gbif_31_butterfly_wgsrpd3_source.py": "7fea5b3849f67c8730cbd246bab67045cb7dff35",
}
PAPER_SOURCE_BLOBS = {
    "aggregate_gbif_global224_wgsrpd3_annual.py": "a06dcb58ffd043623e8227ed4642b17f31c1024e",
    "analyze_gbif_global224_observed_beta.py": "6efdba5c5b47fbb2ef7f1b59876e414004cf76ff",
}
SOURCE_AUDIT = ("provenance/reviewer_defenses/results/"
                "global_colonial_resource_beta_20261031/gbif_global224_annual_v1.audit.json")
GBIF_DOI = "10.15468/dl.umqmsf"
GBIF_KEY = "0013264-261023074918325"
GBIF_ROWS = 14637152
GBIF_ORIGINAL_ARCHIVE_SHA256 = "a1769481e91cb3002c25f69c6d8aa8e6d8d73a012c9495e8394032af02a9df18"


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def imported_relative_modules(script: Path) -> list[str]:
    tree = ast.parse(script.read_text(encoding="utf-8"), filename=str(script))
    results = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            module = node.module.split(".", 1)[0]
            if module.startswith(("analyze_", "audit_", "build_")):
                results.add(module + ".py")
        if isinstance(node, ast.Import):
            for item in node.names:
                module = item.name.split(".", 1)[0]
                if module.startswith(("analyze_", "audit_", "build_")):
                    results.add(module + ".py")
    return sorted(results)


def audit(paper: Path, original: Path) -> dict:
    pinned = {**REQUIRED_SOURCE_BLOBS, **PAPER_SOURCE_BLOBS}
    work = [(name, None) for name in ROOT_SCRIPTS]
    seen = {}
    unresolved = []
    while work:
        name, importing = work.pop()
        if name in seen:
            continue
        f_paper, f_source = paper / "scripts" / name, original / "scripts" / name
        expected = pinned.get(name)
        if f_paper.is_file() and expected == git_blob_sha(f_paper.read_bytes()):
            path, family = f_paper, "paper_original"
        elif f_source.is_file() and expected == git_blob_sha(f_source.read_bytes()):
            path, family = f_source, "source_20261031"
        elif f_paper.is_file() or f_source.is_file():
            unresolved.append({"name": name, "reason": "unrecognized_Git_blob",
                               "imported_from": importing})
            continue
        else:
            unresolved.append({"name": name, "reason": "file_not_in_pinned_checkouts",
                               "imported_from": importing})
            continue
        data = path.read_bytes()
        imports = imported_relative_modules(path)
        seen[name] = {
            "file": name, "origin": family, "git_blob_sha": git_blob_sha(data),
            "sha256": hashlib.sha256(data).hexdigest(), "size": len(data),
            "imports": imports,
        }
        for item in imports:
            work.append((item, name))
    for name, required in REQUIRED_SOURCE_BLOBS.items():
        if name not in seen:
            unresolved.append({"name": name, "reason": "expected_original_helper_not_reached_or_found"})
    for name in ROOT_SCRIPTS:
        if name not in seen:
            unresolved.append({"name": name, "reason": "source_main_entrypoint_missing"})
    source_path = paper / SOURCE_AUDIT
    old_audit = json.loads(source_path.read_text(encoding="utf-8"))
    raw = old_audit["source"]
    if (raw["gbif_doi"] != GBIF_DOI or raw["gbif_download_key"] != GBIF_KEY or
        raw["GBIF_archive_sha256"] != GBIF_ORIGINAL_ARCHIVE_SHA256 or
        old_audit["records"]["rows"] != GBIF_ROWS):
        raise AssertionError("Original GBIF download identity and row count drifted")
    return {
        "schema": "chocho_GEB_v05_gbif224_original_script_dependency_closure_v1",
        "status": ("SOURCE_MODULE_CHAIN_RECOVERED_NOT_RAW_TO_MODEL_RERUN" if not unresolved
                   else "SOURCE_MODULE_CHAIN_INCOMPLETE"),
        "original_paper_snapshot": PAPER_COMMIT_ORIGINAL,
        "original_source_module_commit": SOURCE_COMMIT,
        "shared_source_GBIF_download": {
            "doi": GBIF_DOI, "key": GBIF_KEY,
            "records": GBIF_ROWS, "original_zip_sha256": GBIF_ORIGINAL_ARCHIVE_SHA256,
        },
        "modules": [seen[key] for key in sorted(seen)],
        "missing_source_modules": unresolved,
        "missing_for_full_reproduction": [
            "The original 1.415 GB GBIF occurrence ZIP must be fetched from the same dated download DOI and verified byte-exact before aggregation.",
            "External TDWG level3 shape and the original Lenzner/GloNAF WCVP HOSTS source inputs and original butterfly beta precomputations must be retrieved and verified.",
            "Recomputed annual butterfly records and effort must match source audit (not merely stored fitted scores).",
            "The 155-species RF model train/evaluation and both 224-species primary and geographical-transfer model fits must run end to end and agree with original numeric receipts.",
            "No raw occurrence IDs, observers, coordinates or restricted third-party inputs may be committed in a public Git tree.",
            "An anonymous, licence-cleared, stable reviewer-data route does not yet exist.",
        ],
        "interpretation": "Code lineage only. Independent consumer responses are GBIF recorded communities, not occupancy or causal host use."
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--paper-root", type=Path)
    p.add_argument("--original-root", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        assert git_blob_sha(b"test\n") == hashlib.sha1(b"blob 5\0test\n").hexdigest()
        t=ast.parse("from analyze_a import z\nimport audit_b\nfrom json import loads\n")
        m=set()
        for n in ast.walk(t):
            if isinstance(n, ast.ImportFrom) and n.module and n.module.startswith("analyze_"):
                m.add(n.module + ".py")
            elif isinstance(n, ast.Import):
                m.update(q.name + ".py" for q in n.names if q.name.startswith("audit_"))
        assert m == {"analyze_a.py", "audit_b.py"}
        print("SOURCE_MODULE_AST_SELFTEST_PASS")
        return
    if not a.paper_root or not a.original_root or not a.output:
        p.error("--paper-root --original-root --output required")
    report = audit(a.paper_root, a.original_root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"],
                      "modules": len(report["modules"]),
                      "missing_modules": report["missing_source_modules"]}, ensure_ascii=False))
    if report["missing_source_modules"]:
        raise SystemExit("PENDING_ORIGINAL_CODE_DEPENDENCY")


if __name__ == "__main__":
    main()

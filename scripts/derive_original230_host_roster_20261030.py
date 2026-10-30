#!/usr/bin/env python3
"""Recover original HOSTS catalogue edges for all230 conservatively bridged taxa.

The pinned author ZIP stays in the GitHub runner. Read and verify its members
without extraction, reconcile all239 source counts, and export a bounded derived
roster. No BCE supplementation, occurrence request, plant-presence inference,
raw provider-row export, or archived source-code execution is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "chocho_original230_source_host_roster_20261030_v1"
COMPLETE = "COMPLETE_ORIGINAL230_SOURCE_CATALOGUE_NOT_LOCAL_HOST_USE"
INCOMPLETE = "INCOMPLETE_SOURCE_OR_COHORT_VALIDATION_FAILED"
REPOSITORY = "zuizui0223/chocho"
REPOSITORY_ID = 1390618992
SOURCE_RUN = 39006554609
SOURCE_HEAD = "d53e04fb361701265139ba77b8bca56689da6532"
SOURCE_ARTIFACT = 12090584135
SOURCE_ARTIFACT_NAME = "butterfly-author-anonymous-offline-source-replication-v01"
OUTER_SHA = "a784770f5a215bffb32393c127b66e7320e8ed7ac9b36f921bb9c056f3497a36"
OUTER_SIZE = 3250938
INNER_NAME = "butterfly_resource_novelty_anonymous_OFFLINE_source_v01.zip"
INNER_SHA = "408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
INNER_SIZE = 3255228
HOST_MEMBER = "inputs/original/insect_host_accepted.csv"
METRIC_MEMBER = "inputs/anthropogenic_species_metrics.csv"
MEMBER_PINS = {
    HOST_MEMBER: "0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9",
    METRIC_MEMBER: "7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30",
}
UPSTREAM_PINS = {
    "HOSTS": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
    "WCVP": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
    "TDWG": "52da7828aba9d461dd133c27b3bd7a4407161f54",
}
PLAN = "provenance/archive/host_hindcast/ISLAND_EXACT_NULL_HOST_ROSTER_PLAN_20261030.md"
PLAN_SHA = "ce04b052a37ad4f07d2c81e87916dbf3addac3083a7117ebecea47545a4f99e9"
BRIDGE = "data/exploratory/island_taxonomy_20261029/original239_to_colxr_bridge.csv"
BRIDGE_SHA = "cd8b2cb8e0279a5b940e46267c67f0bb9843d2e1c54f38bc405ec30a95fe6eb2"
FACETS = "data/exploratory/island_expansion_20261029/expanded_island_panel_species_key_facets.csv"
FACETS_SHA = "061f9ce3cf21e4d64572456e0ffbbb88022fbccfe9aa924d5f02062027c97db6"
METRICS = "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
BASELINE_NAME_SHA = "5652125feea89ce2a9228cc19daa6d1e3128784903fa360d15f4a659e4501e1f"
EXPECTED_COUNTS = {
    "original_taxa": 239, "eligible_taxa": 230, "excluded_taxa": 9,
    "original_edges": 2596, "eligible_edges": 2524,
    "baseline_taxa": 93, "baseline_edges": 1440,
}
MANIFEST_NAME = "original230_host_roster_manifest.json"
PLANT_FIELDS = (
    "wcvp_accepted_id", "accepted_name", "accepted_names_json", "family", "families_json",
    "source_name_gate", "ambiguity_flags", "query_eligible", "query_name",
)
BUTTERFLY_FIELDS = (
    "original_name", "current_accepted_key", "current_accepted_name", "current_usage_status",
)
TABLE_FIELDS = {
    "original230_host_edges.csv": BUTTERFLY_FIELDS + (
        "butterfly_match_type", "bridge_response_sha256",
    ) + PLANT_FIELDS,
    "original239_host_count_alignment.csv": BUTTERFLY_FIELDS + (
        "conservative_join_eligible", "exclusion_reasons", "expected_host_count",
        "recovered_host_count", "count_matches", "source_row_count",
    ),
    "plant_taxon_query_candidates.csv": PLANT_FIELDS + (
        "linked_butterfly_count", "linked_original_names_json", "baseline_linked_butterfly_count",
    ),
    "baseline93_host_coverage.csv": BUTTERFLY_FIELDS + (
        "expected_host_count", "recovered_host_count", "query_eligible_host_count",
        "unresolved_host_count", "baseline_record_count", "baseline_island_count",
    ),
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def csv_rows(data, required, label):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError(label + " has missing or duplicate CSV headers")
    if not set(required) <= set(reader.fieldnames):
        raise ValueError(label + " required CSV schema differs")
    rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(label + " has malformed CSV rows")
    return rows


def verify_source_metadata(run, artifact, *, now=None):
    """Validate fresh API identities; export a small receipt, not API payloads."""
    now = now or datetime.now(timezone.utc)
    if (run.get("id"), run.get("head_sha"), run.get("status"), run.get("conclusion")) != (
            SOURCE_RUN, SOURCE_HEAD, "completed", "success"):
        raise ValueError("Pinned source workflow run/head/status changed")
    for key in ("repository", "head_repository"):
        repo = run.get(key, {})
        if (repo.get("id"), repo.get("full_name")) != (REPOSITORY_ID, REPOSITORY):
            raise ValueError("Pinned source repository identity changed")
    if (artifact.get("id"), artifact.get("name"), artifact.get("size_in_bytes"),
            artifact.get("digest"), artifact.get("expired")) != (
            SOURCE_ARTIFACT, SOURCE_ARTIFACT_NAME, OUTER_SIZE, "sha256:" + OUTER_SHA, False):
        raise ValueError("Pinned source artifact identity/digest/size/expiry changed")
    source = artifact.get("workflow_run", {})
    if (source.get("id"), source.get("head_sha"), source.get("repository_id"),
            source.get("head_repository_id")) != (
            SOURCE_RUN, SOURCE_HEAD, REPOSITORY_ID, REPOSITORY_ID):
        raise ValueError("Source artifact workflow provenance changed")
    expiry = datetime.fromisoformat(artifact.get("expires_at", "").replace("Z", "+00:00"))
    if expiry.tzinfo is None or expiry <= now:
        raise ValueError("Pinned source artifact has expired")
    return {
        "repository": REPOSITORY, "run_id": SOURCE_RUN, "run_head_sha": SOURCE_HEAD,
        "run_url": f"https://github.com/{REPOSITORY}/actions/runs/{SOURCE_RUN}",
        "artifact_id": SOURCE_ARTIFACT, "artifact_name": SOURCE_ARTIFACT_NAME,
        "outer_artifact_sha256_from_api_metadata": OUTER_SHA,
        "outer_artifact_bytes_from_api_metadata": OUTER_SIZE,
        "outer_archive_bytes_rehashed_by_this_script": False,
        "artifact_created_utc": artifact.get("created_at"), "artifact_expires_utc": artifact["expires_at"],
        "metadata_verified_utc": now.isoformat(),
    }


def read_verified_archive(path, *, expected_sha256=INNER_SHA, expected_size=INNER_SIZE,
                          expected_members=None):
    """Verify every archive member, then read only original hosts and metrics."""
    expected_members = MEMBER_PINS if expected_members is None else expected_members
    data = Path(path).read_bytes()
    if len(data) != expected_size or sha(data) != expected_sha256:
        raise ValueError("Original inner ZIP size or SHA256 changed")
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        infos = archive.infolist()
        names = [i.filename for i in infos]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate source ZIP member")
        for info in infos:
            name = info.filename
            p = PurePosixPath(name)
            if (not name or "\\" in name or p.is_absolute() or ".." in p.parts
                    or str(p) != name or info.is_dir()
                    or stat.S_ISLNK(info.external_attr >> 16)):
                raise ValueError("Unsafe source ZIP member path or symlink")
        if "SOURCE_MANIFEST.json" not in names or archive.testzip() is not None:
            raise ValueError("Source ZIP manifest missing or CRC check failed")
        manifest_bytes = archive.read("SOURCE_MANIFEST.json")
        manifest = json.loads(manifest_bytes)
        if (manifest.get("schema") != "butterfly_anonymous_offline_resource_package_v0.1"
                or manifest.get("status") != "AUTHOR_REVIEW_ANONYMISED_SOURCE_CONDITIONAL_NOT_A_PUBLIC_DATA_DEPOSIT"
                or manifest.get("upstream_source_pins") != UPSTREAM_PINS
                or manifest.get("sources_conditioned_on_host_catalogue") is not True
                or manifest.get("source_licence_approval_required_before_distribution") is not True):
            raise ValueError("Source manifest schema, scope or upstream pins changed")
        hashes = manifest.get("sha256", {})
        if not isinstance(hashes, dict) or not hashes or set(names) != set(hashes) | {"SOURCE_MANIFEST.json"}:
            raise ValueError("Source manifest does not cover every ZIP member exactly")
        for name, digest in hashes.items():
            if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
                raise ValueError("Malformed source member digest")
            if sha(archive.read(name)) != digest:
                raise ValueError("Archived source member checksum differs")
        if any(hashes.get(name) != digest for name, digest in expected_members.items()):
            raise ValueError("Original host table or source metrics pin changed")
        host_rows = csv_rows(archive.read(HOST_MEMBER),
                            ("insect_species", "accepted_plant_name_id", "accepted_name", "family"),
                            "Original accepted host table")
        metric_rows = csv_rows(archive.read(METRIC_MEMBER),
                              ("species", "resolved_host_species"), "Original239 metrics")
    return {
        "host_rows": host_rows, "metric_rows": metric_rows,
        "provenance": {"inner_zip_sha256": expected_sha256, "inner_zip_bytes": expected_size,
                       "source_manifest_sha256": sha(manifest_bytes),
                       "manifest_members_verified": len(hashes),
                       "input_member_sha256": {name: hashes[name] for name in expected_members},
                       "upstream_source_pins": dict(UPSTREAM_PINS),
                       "archive_members_extracted": 0, "archived_code_executed": False},
    }


def boolean(value, label):
    if value not in ("True", "False"):
        raise ValueError("Invalid Boolean field: " + label)
    return value == "True"


def positive_integer(value, label, *, allow_zero=False):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value):
        raise ValueError("Invalid integer field: " + label)
    number = int(value)
    if number == 0 and not allow_zero:
        raise ValueError("Missing positive value: " + label)
    return number


def butterfly_fields(row):
    return {"original_name": row["original_name"],
            "current_accepted_key": row["current_accepted_key"],
            "current_accepted_name": row["current_accepted_canonical_name"],
            "current_usage_status": row["current_usage_status"]}


def derive_roster(host_rows, metric_rows, bridge_rows, facet_rows, *, expected_counts):
    """Exact source names/IDs determine edges; records only annotate early coverage."""
    metrics = {}
    for row in metric_rows:
        name = row["species"].strip()
        if not name or name in metrics:
            raise ValueError("Duplicate or missing original metric taxon")
        metrics[name] = positive_integer(row["resolved_host_species"], "resolved_host_species")
    bridge = {r["original_name"]: r for r in bridge_rows}
    if len(bridge) != len(bridge_rows) or set(bridge) != set(metrics):
        raise ValueError("Original source239/bridge name identity differs")
    keys = Counter(r["current_accepted_key"] for r in bridge_rows if r["current_accepted_key"])
    eligible, exclusions = {}, {}
    for name, row in bridge.items():
        passed = boolean(row["conservative_join_eligible"], "conservative_join_eligible")
        collision = keys[row["current_accepted_key"]] > 1
        reported_collision = boolean(row["current_accepted_key_collision"], "current_accepted_key_collision")
        if passed and (collision or reported_collision):
            raise ValueError("Eligible butterfly accepted-key collision")
        reasons = set(filter(None, row["rejection_reasons"].split("|")))
        if collision or reported_collision:
            reasons.add("current_accepted_key_collision")
        if passed:
            if (row["match_type"] != "EXACT" or row["current_usage_canonical_name"] != name
                    or row["current_usage_rank"] != "SPECIES" or row["current_accepted_rank"] != "SPECIES"
                    or row["current_usage_status"] not in ("ACCEPTED", "SYNONYM")
                    or not row["current_accepted_key"] or not row["current_accepted_canonical_name"]
                    or reasons or not boolean(row["Papilionoidea_verified"], "Papilionoidea_verified")
                    or not boolean(row["provisional_eligible"], "provisional_eligible")):
                raise ValueError("Conservative butterfly identity gate disagrees with bridge")
            if not re.fullmatch(r"[0-9a-f]{64}", row["response_sha256"]):
                raise ValueError("Butterfly taxonomy response hash missing")
            eligible[name] = row
        else:
            if not reasons:
                raise ValueError("Excluded butterfly concept lacks a documented reason")
            exclusions[name] = "|".join(sorted(reasons))

    by_species = {name: set() for name in metrics}
    source_row_counts = Counter()
    plant_names, plant_families, name_ids = defaultdict(set), defaultdict(set), defaultdict(set)
    missing_names, missing_families = set(), set()
    for row in host_rows:
        name = row["insect_species"].strip()
        pid = row["accepted_plant_name_id"].strip()
        if not re.fullmatch(r"[0-9]+", pid):
            if name in metrics:
                raise ValueError("Missing or non-atomic accepted plant ID for original butterfly")
            continue
        botanical = row["accepted_name"].strip()
        family = row["family"].strip()
        if botanical:
            plant_names[pid].add(botanical)
            name_ids[botanical].add(pid)
        else:
            missing_names.add(pid)
        if family:
            plant_families[pid].add(family)
        else:
            missing_families.add(pid)
        if name in metrics:
            by_species[name].add(pid)
            source_row_counts[name] += 1
    if any(len(by_species[name]) != count for name, count in metrics.items()):
        raise ValueError("Original per-species accepted-host count does not reconcile")

    early_records, early_islands = Counter(), defaultdict(set)
    key_to_name = {r["current_accepted_key"]: name for name, r in eligible.items()}
    seen_early = set()
    for row in facet_rows:
        if row["window"] != "2010_2017":
            continue
        identity = (row["island_key"], row["facet_key"])
        if identity in seen_early:
            raise ValueError("Repeated baseline island/species facet")
        seen_early.add(identity)
        count = positive_integer(row["record_count"], "baseline record_count", allow_zero=True)
        name = key_to_name.get(row["facet_key"])
        if name is not None and count:
            early_records[name] += count
            early_islands[name].add(row["island_key"])
    baseline = set(early_records)
    observed_counts = {
        "original_taxa": len(metrics), "eligible_taxa": len(eligible), "excluded_taxa": len(exclusions),
        "original_edges": sum(map(len, by_species.values())),
        "eligible_edges": sum(len(by_species[n]) for n in eligible),
        "baseline_taxa": len(baseline), "baseline_edges": sum(len(by_species[n]) for n in baseline),
    }
    if observed_counts != expected_counts:
        raise ValueError("Fixed original/eligible/baseline cohort or edge-count totals changed")

    selected_ids = set().union(*(by_species[n] for n in eligible))
    plant_identity = {}
    for pid in sorted(selected_ids):
        names, families = sorted(plant_names[pid]), sorted(plant_families[pid])
        flags = []
        if not names or pid in missing_names:
            flags.append("MISSING_ACCEPTED_NAME")
        if len(names) > 1:
            flags.append("MULTIPLE_ACCEPTED_NAMES")
        if any(not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+", n) for n in names):
            flags.append("NOT_EXACT_SPECIES_BINOMIAL")
        if any(len(name_ids[n]) > 1 for n in names):
            flags.append("NAME_MAPS_TO_MULTIPLE_ACCEPTED_IDS")
        if not families or pid in missing_families:
            flags.append("MISSING_FAMILY")
        if len(families) > 1:
            flags.append("MULTIPLE_FAMILIES")
        query_eligible = not flags and len(names) == 1
        plant_identity[pid] = {
            "wcvp_accepted_id": pid, "accepted_name": names[0] if len(names) == 1 else "",
            "accepted_names_json": json.dumps(names, ensure_ascii=False, separators=(",", ":")),
            "family": families[0] if len(families) == 1 else "",
            "families_json": json.dumps(families, ensure_ascii=False, separators=(",", ":")),
            "source_name_gate": ("EXACT_SOURCE_BINOMIAL_NEEDS_GBIF_TAXON_RESOLUTION" if query_eligible
                                 else "UNRESOLVED_SOURCE_BOTANICAL_IDENTITY"),
            "ambiguity_flags": "|".join(flags), "query_eligible": query_eligible,
            "query_name": names[0] if query_eligible else "",
        }
    edges, alignment, candidates, coverage = [], [], [], []
    linked = defaultdict(list)
    for name in sorted(metrics):
        row = bridge[name]
        alignment.append({**butterfly_fields(row), "conservative_join_eligible": name in eligible,
                          "exclusion_reasons": exclusions.get(name, ""),
                          "expected_host_count": metrics[name], "recovered_host_count": len(by_species[name]),
                          "count_matches": True, "source_row_count": source_row_counts[name]})
        if name not in eligible:
            continue
        for pid in sorted(by_species[name]):
            linked[pid].append(name)
            edges.append({**butterfly_fields(row), "butterfly_match_type": row["match_type"],
                          "bridge_response_sha256": row["response_sha256"], **plant_identity[pid]})
        if name in baseline:
            ready = sum(plant_identity[p]["query_eligible"] for p in by_species[name])
            coverage.append({**butterfly_fields(row), "expected_host_count": metrics[name],
                             "recovered_host_count": len(by_species[name]), "query_eligible_host_count": ready,
                             "unresolved_host_count": len(by_species[name]) - ready,
                             "baseline_record_count": early_records[name],
                             "baseline_island_count": len(early_islands[name])})
    for pid in sorted(selected_ids):
        candidates.append({**plant_identity[pid], "linked_butterfly_count": len(linked[pid]),
                           "linked_original_names_json": json.dumps(linked[pid], ensure_ascii=False,
                                                                    separators=(",", ":")),
                           "baseline_linked_butterfly_count": len(baseline & set(linked[pid]))})
    flag_counts = Counter(flag for r in candidates for flag in r["ambiguity_flags"].split("|") if flag)
    return {
        "tables": {"original230_host_edges.csv": edges, "original239_host_count_alignment.csv": alignment,
                   "plant_taxon_query_candidates.csv": candidates, "baseline93_host_coverage.csv": coverage},
        "summary": {"original_taxa": len(metrics), "eligible_taxa": len(eligible),
                    "excluded_taxa": len(exclusions), "original_catalog_edges": observed_counts["original_edges"],
                    "exported_catalog_edges": len(edges), "unique_source_plant_ids": len(candidates),
                    "candidate_exact_name_query_ids": sum(r["query_eligible"] for r in candidates),
                    "unresolved_candidate_ids": sum(not r["query_eligible"] for r in candidates),
                    "baseline_taxa": len(baseline), "baseline_catalog_edges": observed_counts["baseline_edges"],
                    "baseline_taxa_with_all_catalog_edges_recovered": len(coverage),
                    "baseline_taxa_with_unresolved_plant_names": sum(r["unresolved_host_count"] > 0 for r in coverage),
                    "baseline_original_names_sha256": sha(("\n".join(sorted(baseline)) + "\n").encode())},
        "taxonomy_audit": {"butterfly_current_usage_statuses": dict(sorted(Counter(
                               r["current_usage_status"] for r in eligible.values()).items())),
                           "excluded_original_concepts": [{"original_name": n, "reasons": exclusions[n]}
                                                          for n in sorted(exclusions)],
                           "plant_ambiguity_flag_counts": dict(sorted(flag_counts.items())),
                           "plant_identity_conflicts_checked_across_original_source_table": True},
    }


def export_roster(result, output, *, provenance):
    """Exact field and filename allowlists gate every output before any writing."""
    output = Path(output)
    tables = result["tables"]
    if set(tables) != set(TABLE_FIELDS):
        raise ValueError("Derived output table allowlist differs")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Derived output directory must be empty")
    encoded = {}
    for name, fields in TABLE_FIELDS.items():
        rows = tables[name]
        if any(set(row) != set(fields) for row in rows):
            raise ValueError("Derived output column allowlist differs")
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        encoded[name] = stream.getvalue().encode("utf-8")
    manifest = {
        "schema": SCHEMA, "status": COMPLETE, "executed_utc": utc_now(),
        "analysis_date_JST": "2026-10-30", "summary": result["summary"],
        "taxonomy_audit": result["taxonomy_audit"], "provenance": provenance,
        "output_sha256": {n: sha(data) for n, data in sorted(encoded.items())},
        "output_allowlist": sorted(set(TABLE_FIELDS) | {MANIFEST_NAME}),
        "catalogue_scope": "Original source HOSTS accepted IDs only; no BCE addition; all230 conservative taxa.",
        "baseline_scope": "Only positive 2010_2017 expanded96 species facets define the93 coverage taxa; no late selection.",
        "query_eligibility_definition": "Unambiguous exact source binomial and family; separate GBIF taxon resolution is still required.",
        "no_occurrence_queries_performed": True, "no_raw_provider_rows_or_coordinates_exported": True,
        "plant_occurrence_taxonomy_resolution_performed": False,
        "interpretation": "Source catalogue associations and identity coverage only. These are not island plant presence, dated introduction, actual feeding, recruitment, or occurrence absence. Current butterfly synonym mappings do not independently validate unchanged host concepts.",
    }
    encoded[MANIFEST_NAME] = (json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()
    output.mkdir(parents=True, exist_ok=True)
    for name, data in encoded.items():
        (output / name).write_bytes(data)
    return manifest


def checked_local(root, relative, expected):
    data = (root / relative).read_bytes()
    if sha(data) != expected:
        raise ValueError("Pinned local input hash differs: " + relative)
    return data


def execution_receipt():
    keys = ("GITHUB_SHA", "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "GITHUB_JOB",
            "GITHUB_WORKFLOW", "GITHUB_WORKFLOW_REF", "GITHUB_REF", "GITHUB_REPOSITORY")
    return {key: os.environ[key] for key in keys if os.environ.get(key)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--run-metadata", type=Path, required=True)
    parser.add_argument("--artifact-metadata", type=Path, required=True)
    parser.add_argument("--verify-source-metadata", action="store_true")
    parser.add_argument("--source-archive-directory", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    stage = "plan_and_source_metadata"
    try:
        checked_local(args.repo_root, PLAN, PLAN_SHA)
        source_metadata = verify_source_metadata(json.loads(args.run_metadata.read_text()),
                                                 json.loads(args.artifact_metadata.read_text()))
        if args.verify_source_metadata:
            print("ORIGINAL230_SOURCE_METADATA_VERIFIED=" + json.dumps(source_metadata, sort_keys=True))
            return 0
        if args.source_archive_directory is None or args.output_dir is None:
            raise ValueError("Archive directory and output directory are required for derivation")
        stage = "verified_original_archive"
        archives = list(args.source_archive_directory.rglob(INNER_NAME))
        if len(archives) != 1 or not archives[0].is_file() or archives[0].is_symlink():
            raise ValueError("Exact source inner ZIP is missing or duplicated")
        source = read_verified_archive(archives[0])
        stage = "pinned_local_cohorts"
        metric_bytes = checked_local(args.repo_root, METRICS, MEMBER_PINS[METRIC_MEMBER])
        local_metrics = csv_rows(metric_bytes, ("species", "resolved_host_species"), "Local original239 metrics")
        if local_metrics != source["metric_rows"]:
            raise ValueError("Local and archive original239 metric tables differ")
        bridge = csv_rows(checked_local(args.repo_root, BRIDGE, BRIDGE_SHA),
                          ("original_name", "conservative_join_eligible"), "Conservative239 bridge")
        facets = csv_rows(checked_local(args.repo_root, FACETS, FACETS_SHA),
                          ("island_key", "window", "facet_key", "record_count"), "Expanded96 facets")
        stage = "all239_reconciliation_and_all230_derivation"
        result = derive_roster(source["host_rows"], source["metric_rows"], bridge, facets,
                               expected_counts=EXPECTED_COUNTS)
        if result["summary"]["baseline_original_names_sha256"] != BASELINE_NAME_SHA:
            raise ValueError("Baseline93 exact original-name roster differs")
        provenance = {"source_artifact": source_metadata, "source_archive": source["provenance"],
                      "local_input_sha256": {PLAN: PLAN_SHA, BRIDGE: BRIDGE_SHA, FACETS: FACETS_SHA,
                                             METRICS: MEMBER_PINS[METRIC_MEMBER]},
                      "script_sha256": sha(Path(__file__).read_bytes()),
                      "execution": execution_receipt(), "python_version": sys.version.split()[0]}
        stage = "derived_only_export"
        manifest = export_roster(result, args.output_dir, provenance=provenance)
        print("ORIGINAL230_ROSTER_MANIFEST_SUMMARY=" + json.dumps(
            {k: manifest[k] for k in ("schema", "status", "executed_utc", "summary", "provenance", "output_sha256")},
            ensure_ascii=False, sort_keys=True, separators=(",", ":")), flush=True)
        return 0
    except (ValueError, OSError, KeyError, TypeError, zipfile.BadZipFile, csv.Error) as exc:
        failure = {"schema": SCHEMA, "status": INCOMPLETE, "failed_stage": stage,
                   "executed_utc": utc_now(), "error_type": type(exc).__name__,
                   "message": str(exc)[:300], "qualified_host_roster": False,
                   "no_occurrence_queries_performed": True, "execution": execution_receipt()}
        if args.output_dir is not None and (not args.output_dir.exists() or (
                args.output_dir.is_dir() and not any(args.output_dir.iterdir()))):
            args.output_dir.mkdir(parents=True, exist_ok=True)
            (args.output_dir / MANIFEST_NAME).write_text(json.dumps(failure, indent=2) + "\n")
        print("ORIGINAL230_ROSTER_INCOMPLETE=" + json.dumps(failure, sort_keys=True), flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

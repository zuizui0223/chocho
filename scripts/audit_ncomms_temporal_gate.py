#!/usr/bin/env python3
"""Source-locked Nature Communications feasibility screen for temporal host × access.

This script audits *archived* 2001–2010 host-record support for the later
(2011–2022) butterfly record-index cohort, with original 225-species identities.
It does NOT fit an ecological model or mistake non-records for true absences.
All data inputs are already tracked in the independent chocho source branch.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = "analysis/ncomms_temporal_access_20261108"
SOURCE = {
    "species": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_all225.csv",
        "4f3ebbdd5f540ae0eeedf5697a1922f166de7076",
    ),
    "dated_validation": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_validation.json",
        "1109c17b4c0de938e6b6defe20d05db3c5b2a88b",
    ),
    "dated_family_null": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_host_family_null_dated_audit.json",
        "24ba0c8a66a788ec9a463311f95c0b0b9d3388ce",
    ),
    "early_inventory": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_audit.json",
        "95a8ea7dc2096611837ebbd66aa27025b7fc9c96",
    ),
    "early_island_roster": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_islands.csv",
        "143f9956ce4f871c7a97b4777ab69d416963c2bb",
    ),
    "early_general_island": (
        "results/island51_pre2011_public_general_flora_20261107/source/GBIF_original51_early_general_flora_islands.csv",
        "048fa6d644f5579879620deb394ce74ace653e26",
    ),
    "early_general_ports": (
        "results/island51_pre2011_public_general_flora_20261107/source/GBIF_original194_port_2001_2010_coverage.csv",
        "32217fc3f8d0773c679cd85a3420d890e6967200",
    ),
    "dated_comparison": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_comparison.csv",
        "15b2282b0cd255a2f3d518e495727353bee8b1e0",
    ),
}


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def read_pinned(root: Path = ROOT) -> tuple[dict, dict]:
    records, provenance = {}, {}
    for key, (relative, fixed_sha) in SOURCE.items():
        source = root / relative
        raw = source.read_bytes()
        measured_sha = git_blob_sha(raw)
        if measured_sha != fixed_sha:
            raise ValueError(f"Original source blob does not match pre-locked receipt: {relative}")
        if relative.endswith(".csv"):
            data = list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
        else:
            data = json.loads(raw.decode("utf-8"))
        records[key] = data
        provenance[key] = {"repo_path": relative, "source_git_blob_sha1": measured_sha,
                           "raw_sha256": hashlib.sha256(raw).hexdigest()}
    return records, provenance


def fixed_integer(raw: str | float | int, *, label: str) -> int:
    if raw is None or raw == "":
        raise ValueError(f"Missing integer {label}")
    value = float(raw)
    if not math.isfinite(value) or value < 0 or not value.is_integer():
        raise ValueError(f"Invalid integer {label}: {raw!r}")
    return int(value)


def summarize_species(rows: list[dict]) -> dict:
    if len(rows) != 225:
        raise ValueError(f"Original candidate universe changed: {len(rows)} != 225")
    seen_names, seen_keys = set(), set()
    ev, observed_pairs, hits, matched_species = 0, 0, 0, []
    for row in rows:
        species = row["species"].strip()
        key = row["taxon_key"].strip()
        if not species or not key or species in seen_names or key in seen_keys:
            raise ValueError("Empty/duplicate original butterfly name or GBIF taxonKey")
        seen_names.add(species)
        seen_keys.add(key)
        eligible = row["comparison_evaluable"].strip().lower() == "true"
        if row["comparison_evaluable"].strip().lower() not in {"true", "false"}:
            raise ValueError(f"Invalid evaluation flag for {species}")
        if not eligible:
            if any(row[k].strip() for k in
                   ("both_new_pairs", "matching_early_pairs", "fraction_matching_early")):
                raise ValueError(f"Excluded species has non-missing counted pairs: {species}")
            continue
        ev += 1
        n = fixed_integer(row["both_new_pairs"], label=species + " later records")
        m = fixed_integer(row["matching_early_pairs"], label=species + " previous-host matches")
        if n == 0 or m > n:
            raise ValueError(f"Invalid temporal-support denominator for {species}")
        fraction = float(row["fraction_matching_early"])
        if not math.isfinite(fraction) or not math.isclose(fraction, m/n, abs_tol=1e-12):
            raise ValueError(f"Incorrect dated-host match fraction for {species}")
        observed_pairs += n
        hits += m
        if m:
            matched_species.append({"butterfly": species,
                                    "later_record_pairs": n, "with_prior_host_record": m})
    return {"original_candidate_butterflies": len(rows),
            "butterflies_with_later_source_pairs": ev,
            "later_indexed_butterfly_species_island_pairs": observed_pairs,
            "later_pairs_with_recorded_earlier_matching_host": hits,
            "butterfly_species_with_verified_prior_host_overlap": len(matched_species),
            "exact_butterflies_with_earlier_recorded_host_overlap": matched_species,
            "fraction_of_later_record_pairs_with_prior_host_record": hits/observed_pairs,
            "later_pair_records_without_observed_prior_host_match": observed_pairs-hits}


def audit(records: dict, provenance: dict) -> dict:
    species = summarize_species(records["species"])
    early, control, null = (records[x] for x in
                            ("early_inventory", "dated_validation", "dated_family_null"))
    n, m = (species[x] for x in
            ("later_indexed_butterfly_species_island_pairs",
             "later_pairs_with_recorded_earlier_matching_host"))
    if control["status"] != "PASS_DATED_MATCHES_RECOMPUTED":
        raise ValueError("Dated-host validation is not PASS")
    for key, val in [("records_checked", n), ("matched_pairs", m),
                     ("matched_species", species["butterfly_species_with_verified_prior_host_overlap"])]:
        if int(control[key]) != val:
            raise ValueError(f"Disagreement with independent dated-host validation: {key}")
    if not control["all225_accounted"] or int(control["changed_count_ports"]) != 0:
        raise ValueError("Original 225 candidate population, or original port retrieval, changed")
    for key, expected in [("planned_ports", 194), ("complete_ports", 194),
                          ("completely_resolved_islands", 51), ("host_port_pairs", 7967)]:
        if int(early[key]) != expected:
            raise ValueError(f"Historical botanical data collection changed: {key}")
    if int(early["live_query_completed_ports"]) + int(early["cached_zero_ports"]) != 194:
        raise ValueError("Botanical inventory source retrieval fails coverage total")

    if int(null["candidate_species"]) != 225 or int(null["both_new_species"]) != species["butterflies_with_later_source_pairs"]:
        raise ValueError("Original family-preserving null candidate/species counts changed")
    if int(null["both_new_species_pairs"]) != n or int(null["draws"]) != 199:
        raise ValueError("Original degree-preserving dated null source changed")
    if not math.isclose(float(null["observed_fraction"]), m/n, abs_tol=1e-13):
        raise ValueError("Archived null observed denominator disagrees with dated-host table")

    by_name = {row["plant_window"]: row for row in records["dated_comparison"]}
    if set(by_name) != {"Plants2001-2024", "Plants2001-2010"}:
        raise ValueError("Expected two genuinely different temporal windows")
    modern, previous = (by_name[k] for k in ("Plants2001-2024", "Plants2001-2010"))
    if not math.isclose(float(previous["observed"]), m/n, abs_tol=1e-13):
        raise ValueError("Pre-2011 source comparison changed")
    n_pooled = round(float(modern["observed"])*n)
    if n_pooled != 208:
        raise ValueError("Pooled plant record-coverage comparison changed")
    if int(previous["draws_at_least_observed"]) != int(null["draws_at_least_observed"]) or (
       int(previous["draws_at_least_species_equal"]) != int(null["draws_at_least_observed_species_equal"])):
        raise ValueError("Degree/family null output disagrees with shared comparison file")

    source_islands = records["early_island_roster"]
    historical = records["early_general_island"]
    port_rows = records["early_general_ports"]
    if len(source_islands) != 51 or len(historical) != 51 or len(port_rows) != 194:
        raise ValueError("Changed original 51-island / 194-port source population")
    ids = {row["landmass_id"] for row in source_islands}
    if len(ids) != 51 or ids != {row["island_id"] for row in historical}:
        raise ValueError("Historical botanical snapshots have different island IDs")
    if sum(fixed_integer(z["planned_ports"], label="roster ports") for z in source_islands) != 194:
        raise ValueError("Original island port counts no longer sum to 194")
    if any(z["all_ports_complete"] != "True" or z["planned_ports"] != z["complete_ports"]
           for z in source_islands):
        raise ValueError("Some original source island's botanical ports are incomplete")
    if {z["port_id"] for z in port_rows}.__len__() != 194:
        raise ValueError("Duplicate original port in new public general flora")
    observed_plant_ports = sum(int(z["live_GBIF_2001_2010_plant_records"]) > 0 for z in port_rows)
    indexed_plant_islands = sum(int(z["unique_early_indexed_general_plant_speciesKeys"]) > 0 for z in historical)
    if observed_plant_ports != 149 or indexed_plant_islands != 45:
        raise ValueError("Public 2001–2010 **general** flora coverage changed")

    # A block because original HOSTS × port × year identities are not available
    # in the reproducible Git tree. The 194-port general flora is not a substitute.
    # This is a conservative support screen, not a power test.
    heuristic = {"prior_record_supported_later_cases": 50,
                 "different_butterfly_species_with_prior_matches": 10,
                 "different_island_groups_with_prior_matches": 10}
    later_case_support_passes = (
        m >= heuristic["prior_record_supported_later_cases"] and
        len(species["exact_butterflies_with_earlier_recorded_host_overlap"]) >=
        heuristic["different_butterfly_species_with_prior_matches"]
    )
    # There are at most five distinct islands with matched prior-host cases.
    # Do not invent which of the 51 they are: the public aggregated table
    # does not contain the island IDs of the five matches.
    island_criterion_possible = m >= heuristic["different_island_groups_with_prior_matches"]
    readiness = ("MAY_PROCEED_TO_DETECTION_AWARE_DESIGN_AFTER_RAW_SOURCE_RECOVERY"
                 if later_case_support_passes and island_criterion_possible else
                 "BLOCKED_PRIOR_HOST_SUPPORT_AND_MISSING_AT_RISK_SOURCE")

    return {
        "schema": "chocho_NC_temporal_resource_access_source_audit_v1",
        "status": "ARCHIVED_SOURCE_COUNTS_RECOMPUTED_NO_ECOLOGICAL_MODEL",
        "publication_target": "Nature Communications (target only, not submitted)",
        "primary_test_readiness": readiness,
        "source_pins": provenance,
        "prior_host_record_support": species,
        "source_coverage": {
            "original_source_landmass_groups": 51,
            "original_source_ports": 194,
            "original_2001_2010_host_record_inventory_completed_ports": int(early["complete_ports"]),
            "archived_original_2001_2010_host_port_pairs": int(early["host_port_pairs"]),
            "historical_general_plant_port_queries_with_indexed_records": observed_plant_ports,
            "historical_general_plant_islands_with_indexed_speciesKeys": indexed_plant_islands,
            "historical_general_plant_islands_without_indexed_speciesKeys": 51-indexed_plant_islands,
            "historical_general_plants_are_not_confirmed_butterfly_hosts": True,
        },
        "retrospective_vs_prior_comparison": {
            "later_record_pairs_with_matching_host_2001_2024": n_pooled,
            "later_record_pairs_with_matching_host_2001_2010": m,
            "later_record_pairs_in_common_denominator": n,
            "extra_matches_in_pooled_window_without_verified_prior_support": n_pooled-m,
            "NOT_A_TEMPORAL_INTRODUCTION_ESTIMATE": True,
        },
        "dated_host_family_preserving_null": {
            "draws": int(null["draws"]),
            "observed_pair_match_fraction": float(null["observed_fraction"]),
            "null_mean_pair_match_fraction": float(null["null_mean_fraction"]),
            "draws_matching_or_exceeding_observed_pair_fraction": int(null["draws_at_least_observed"]),
            "draws_matching_or_exceeding_observed_species_equal": int(null["draws_at_least_observed_species_equal"]),
            "causal_test": False,
        },
        "screening_only_heuristic_not_formal_power": heuristic,
        "screening_passed": bool(later_case_support_passes and island_criterion_possible),
        "archived_summary_cannot_identify_the_five_matching_islands": True,
        "current_source_cannot_demonstrate": [
            "Plant introduction year", "Prior host's reproductive fitness effect",
            "Dispersal via specific ships", "Real butterfly arrival or breeding",
            "Complete at-risk species × island non-event population",
            "Species-specific access × prior-host temporal interaction",
        ],
        "needed_before_model": [
            "Publicly auditable, dated accepted host-by-port and island species identities",
            "Complete pre-baseline and follow-up butterfly species × island risk-set indexing with explicit missingness",
            "Pre-2011 host and butterfly detection-effort controls plus climate/area",
            "Fixed species- and island-block holdout, structural non-host and source-null controls",
            "Independent spatial or taxonomic validation and sensitivity to source indexing changes",
        ],
    }


def render_markdown(report: dict) -> str:
    s=report["prior_host_record_support"]
    a=report["source_coverage"]
    c=report["retrospective_vs_prior_comparison"]
    n=s["later_indexed_butterfly_species_island_pairs"]
    m=s["later_pairs_with_recorded_earlier_matching_host"]
    names=", ".join(z["butterfly"] for z in s["exact_butterflies_with_earlier_recorded_host_overlap"])
    return f"""# Nature Communications temporal support diagnostic — original sources only

**Result: {report["primary_test_readiness"]}**. The source audit passed, but this is **not a new fitted predictive result**, an ecological absence estimate, or a rejection of the existing island/port positives.

- Original panel: **{s["original_candidate_butterflies"]}** butterflies; **{s["butterflies_with_later_source_pairs"]}** with later 2011–2022 butterfly pairs; **{n}** original butterfly×island later-index records.
- Verified *recorded* pre-2011 matching host: **{m}/{n} ({100*m/n:.2f}%)**, across **{s["butterfly_species_with_verified_prior_host_overlap"]}** butterfly species: {names}. Remaining {n-m} are unmatched in this record-based screen, **not known host absences**.
- The same later-pair population with **pooled 2001–2024** host records gives **{c["later_record_pairs_with_matching_host_2001_2024"]}/{n}** observed host correspondences. The **{c["extra_matches_in_pooled_window_without_verified_prior_support"]}** extra pooled-window matches cannot be claimed as time-prior hosts.
- Original old-host-record acquisition: **{a["original_2001_2010_host_record_inventory_completed_ports"]}/{a["original_source_ports"]}** WPI port buffers and all **{a["original_source_landmass_groups"]}** island groups were audited. The separately reconstructed public early **general** plant flora has index evidence in **{a["historical_general_plant_port_queries_with_indexed_records"]}/{a["original_source_ports"]}** ports and **{a["historical_general_plant_islands_with_indexed_speciesKeys"]}/51** island groups; this is not dated butterfly host identity.
- The original family-preserving null found **{report["dated_host_family_preserving_null"]["draws_matching_or_exceeding_observed_pair_fraction"]}/199** observational surrogates at least as high as the dated pair-matched fraction. No pair-independent or causal p-value is inferred.
- The proposed 50-case/10-species/10-island source-support benchmark is an **exploratory screening rule, not a statistical power calculation**, and is **not met**.

## Required next evidence, without changing positive historical findings

{chr(10).join("- "+x for x in report["needed_before_model"])}

## Scope
Preserve the **239-butterfly WGSRPD3** global source separately from this **225-butterfly/51-island** observational source. Retain the previously established resource and port-flora associations and the exploratory marine-circuit predictive increments, but do not convert them into verified plant movement, host use or butterfly colonisation. Nature Communications remains the **target journal**, while model-ready evidence is not yet established.
"""


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir",type=Path,default=ROOT/OUTDIR/"screen_output")
    args=parser.parse_args()
    data, pins=read_pinned()
    result=audit(data,pins)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"temporal_source_gate.json").write_text(json.dumps(result,indent=2)+"\n")
    (args.out_dir/"SUMMARY.md").write_text(render_markdown(result))
    print(json.dumps({"audit":result["status"],
                      "readiness":result["primary_test_readiness"],
                      "butterflies":result["prior_host_record_support"]["original_candidate_butterflies"],
                      "later_pairs":result["prior_host_record_support"]["later_indexed_butterfly_species_island_pairs"],
                      "prior_host_matches":result["prior_host_record_support"]["later_pairs_with_recorded_earlier_matching_host"],
                      "output":str(args.out_dir)},indent=2))


if __name__=="__main__":
    main()

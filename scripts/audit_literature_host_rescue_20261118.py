#!/usr/bin/env python3
"""Audit independent literature larval-host reports against name-reconciled NHM HOSTS.

This is a PURPOSEFULLY TARGETED source-recovery audit on previously HOSTS-unmatched
island-endemic butterflies. The resulting union is NOT an unbiased estimate of
biological host-knowledge coverage or a cross-stratum sampling comparison.
"""
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

PARENT = Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
HERE = Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
FROZEN42 = PARENT / "independent_42_species_vs_frozen239.csv"
RESCUE = HERE / "independent_larval_host_rescue_10.csv"
CANARY_BIOTA = HERE / "canary_BIOTA_source_status_20.csv"
EVIDENCE_TYPES = {
    "secondary_peer_review_conservation",
    "first_hand_field_natural_history",
    "peer_reviewed_literature_synthesis",
    "curated_natural_history",
}
RESCUE_COLUMNS = [
    "independent_checklist_species", "archipelago",
    "larval_host_genera", "source_url", "evidence_type", "evidence_note",
]


def load(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def yes(value):
    if str(value) not in {"0", "1"}:
        raise ValueError("Binary audit flag expected: " + str(value))
    return int(value)


def run(frozen, reconciled, rescue, biota, output):
    checklist = load(frozen)
    hosts = load(reconciled)
    lit = load(rescue)
    island_status = load(biota)
    if len(checklist) != 42 or len(hosts) != 42:
        raise ValueError("Independent 42 butterfly denominator changed")
    if len(lit) != 10 or len(island_status) != 20:
        raise ValueError("Independent rescue/BIOTA evidence sheet count changed")
    by_checklist = {r["scientific_name"]: r for r in checklist}
    by_hosts = {r["scientific_name"]: r for r in hosts}
    biota_names = {r["independent_checklist_name"]: r for r in island_status}
    if (len(by_checklist), len(by_hosts), len(biota_names)) != (42, 42, 20):
        raise ValueError("Duplicate independent butterfly names or origin entries")
    if set(by_hosts) != set(by_checklist):
        raise ValueError("Reconciled names differ from independently listed 42 butterfly species")
    if {r["independent_checklist_species"] for r in lit} & {
        n for n, h in by_hosts.items()
        if yes(h["hosts_genus_after_guarded_synonyms"])
    }:
        raise ValueError("Rescued literature host taxon overlaps name-reconciled HOSTS")
    if len({r["independent_checklist_species"] for r in lit}) != len(lit):
        raise ValueError("Duplicated independent literature taxon")
    if sum(yes(r["hosts_named_genus"]) for r in checklist) != 19:
        raise ValueError("Original exact name baseline unexpected")
    if sum(yes(r["hosts_genus_after_guarded_synonyms"]) for r in hosts) != 21:
        raise ValueError("Taxonomically guarded HOSTS baseline unexpected")
    for name, row in by_hosts.items():
        old = by_checklist[name]
        if old["checklist_family"].lower() != row["family"].lower():
            raise ValueError("Taxonomic family changed for " + name)
        if yes(old["hosts_named_genus"]) != yes(row["hosts_exact_genus"]):
            raise ValueError("HOSTS source exact genus flag changed for " + name)
        if yes(old["in_frozen239"]) != yes(row["in_original239"]):
            raise ValueError("Original GEB239 inclusion changed for " + name)
    by_literature = {}
    for record in lit:
        name = record["independent_checklist_species"]
        if name not in by_checklist:
            raise ValueError("Independent source taxon not in checklist: " + name)
        row = by_checklist[name]
        if yes(row["hosts_named_genus"]) or yes(by_hosts[name]["hosts_genus_after_guarded_synonyms"]):
            raise ValueError("Literature source counted as fresh but already in HOSTS: " + name)
        if record["archipelago"] not in row["archipelagos"].split("; "):
            raise ValueError("Source region not in independent checklist: " + name)
        if record["archipelago"] == "Madeira":
            if row["madeira_origin"] != "island_endemic_species":
                raise ValueError("Madeira rescue does not have independent endemic origin")
        elif record["archipelago"] == "Canary Islands":
            origin = biota_names.get(name, {}).get("canary_origin_stratum")
            if origin != "endemic_confirmed":
                raise ValueError("Canary rescue lacks independent BIOTA endemic origin: " + name)
        else:
            raise ValueError("Literature rescue comes from unverified origin group")
        if record["evidence_type"] not in EVIDENCE_TYPES:
            raise ValueError("Unknown evidence category for " + name)
        if not record["source_url"].startswith("https://"):
            raise ValueError("No durable cited source URL for " + name)
        if not record["evidence_note"].strip():
            raise ValueError("No description to defend host type: " + name)
        genera = record["larval_host_genera"].split(";")
        if not genera or len(genera) != len(set(genera)):
            raise ValueError("No genera or duplicate names for " + name)
        if any(not re.fullmatch(r"[A-Z][a-zA-Z-]*", genus) for genus in genera):
            raise ValueError("Non-binomial or invalid genus in " + name)
        by_literature[name] = record
    per_taxon = []
    for name in sorted(by_checklist):
        old = by_checklist[name]
        new = by_hosts[name]
        lit_rec = by_literature.get(name)
        has_reconciled_host = yes(new["hosts_genus_after_guarded_synonyms"])
        row = {
            "scientific_name": name,
            "family": old["checklist_family"],
            "archipelagos": old["archipelagos"],
            "hosts_exact_named_genus": yes(old["hosts_named_genus"]),
            "hosts_guarded_synonym_named_genus": has_reconciled_host,
            "targeted_outside_HOSTS_larval_host_genus": int(bool(lit_rec)),
            "combined_confirmed_named_host_genus": int(bool(lit_rec) or bool(has_reconciled_host)),
            "in_original_GEB239": yes(old["in_frozen239"]),
            "literature_host_genera": (lit_rec["larval_host_genera"] if lit_rec else ""),
            "literature_source_url": (lit_rec["source_url"] if lit_rec else ""),
            "literature_evidence_type": (lit_rec["evidence_type"] if lit_rec else ""),
            "madeira_origin": old["madeira_origin"],
            "canary_origin": biota_names.get(name, {}).get("canary_origin_stratum", ""),
        }
        per_taxon.append(row)

    def subset_stats(names):
        found = [r for r in per_taxon if r["scientific_name"] in names]
        if len(found) != len(names):
            raise ValueError("Cohort name absent after join")
        return {
            "n": len(found),
            "exact_HOSTS": sum(r["hosts_exact_named_genus"] for r in found),
            "guarded_taxonomy_HOSTS": sum(r["hosts_guarded_synonym_named_genus"] for r in found),
            "targeted_additional_literature": sum(r["targeted_outside_HOSTS_larval_host_genus"] for r in found),
            "combined_confirmed_minimum": sum(r["combined_confirmed_named_host_genus"] for r in found),
            "in_original_GEB239": sum(r["in_original_GEB239"] for r in found),
        }
    madeira = {n for n, r in by_checklist.items() if r["madeira_origin"] == "island_endemic_species"}
    canary_endemic = {n for n, r in biota_names.items() if r["canary_origin_stratum"] == "endemic_confirmed"}
    canary_native = {n for n, r in biota_names.items() if r["canary_origin_stratum"] in {
        "native_nonendemic_secure", "native_nonendemic_probable"
    }}
    azores_endemic = {"Hipparchia azorina", "Hipparchia miguelensis"}
    strata = {
        "Madeira_extant_endemic": subset_stats(madeira),
        "Canary_BIOTA_endemic_13": subset_stats(canary_endemic),
        "Canary_BIOTA_nonendemic_native_7": subset_stats(canary_native),
        "Azores_Hipparchia_endemic_2": subset_stats(azores_endemic),
    }
    if strata["Madeira_extant_endemic"]["n"] != 3 or strata["Madeira_extant_endemic"]["targeted_additional_literature"] != 3:
        raise ValueError("Not all three Madeira endemic literature confirmations preserved")
    if strata["Canary_BIOTA_endemic_13"]["n"] != 13 or strata["Canary_BIOTA_endemic_13"]["targeted_additional_literature"] != 7:
        raise ValueError("Seven independent Canary rescue confirmations missing")
    if strata["Canary_BIOTA_nonendemic_native_7"]["n"] != 7 or strata["Canary_BIOTA_nonendemic_native_7"]["exact_HOSTS"] != 4:
        raise ValueError("Independent Canary native controls changed")
    if strata["Azores_Hipparchia_endemic_2"]["guarded_taxonomy_HOSTS"] != 2:
        raise ValueError("Azores endemic HOSTS control changed")
    combined = sum(r["combined_confirmed_named_host_genus"] for r in per_taxon)
    if combined != 31:
        raise ValueError("Unexpected total source-aware lower-bound union")
    summary = {
        "schema": "chocho_targeted_larval_host_literature_rescue_v1",
        "independent_unique_butterflies": len(per_taxon),
        "original_exact_NHM_HOSTS_named_genus": 19,
        "guarded_synonym_NHM_HOSTS_named_genus": 21,
        "independent_targeted_larval_host_rescues_not_in_HOSTS": 10,
        "combined_ever_documented_named_host_genus_minimum": combined,
        "combined_minimum_percent_of_independent_42": round(100 * combined / 42, 4),
        "targeted_search_only": True,
        "source_observation_types": dict(Counter(x["evidence_type"] for x in lit)),
        "strata": strata,
        "source_files": {
            "independent_42": str(frozen),
            "guarded_synonym_HOSTS_42": str(reconciled),
            "independent_origin_BIOTA_20": str(biota),
            "targeted_endemic_literature_10": str(rescue),
        },
        "inference_boundary": [
            "A lower bound on distinct butterfly species with a named reported larval-host genus in HOSTS and the ten targeted independently referenced sources, NOT an unbiased true host-known rate.",
            "The ten added host-reported endemic butterflies were searched purposefully because they lacked HOSTS records; native and endemic groups were not searched symmetrically.",
            "The source evidence categories range from firsthand wild-larva observations to secondary curated natural-history summaries.",
            "Synonym-corrected HOSTS identity does not resolve conflicting island-endemic species concepts. Names still taxonomically unresolved are retained.",
            "Larval host-genus knowledge alone does not establish native-range plant geography, local resource availability, larval performance, or island colonization.",
            "Prior GEB 239-species known-host opportunity results remain original, unchanged, and selected-population conditional.",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    with (output / "confirmed_named_host_per_butterfly_42.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(per_taxon[0]))
        writer.writeheader()
        writer.writerows(per_taxon)
    (output / "literature_host_rescue_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print("LITERATURE_HOST_RESCUE_SOURCE_AUDIT_PASSED")
    print("LITERATURE_HOST_RESCUE_SUMMARY=" + json.dumps({
        "host_name_reconciled": 21, "targeted_literature_additions": 10,
        "union_confirmed_lower_bound": combined,
        "out_of_independent_42": 42, "cohorts": strata
    }, ensure_ascii=False))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frozen-42", type=Path, default=FROZEN42)
    ap.add_argument("--reconciled-42", type=Path, required=True)
    ap.add_argument("--rescue-10", type=Path, default=RESCUE)
    ap.add_argument("--canary-biota-20", type=Path, default=CANARY_BIOTA)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    run(args.frozen_42, args.reconciled_42, args.rescue_10, args.canary_biota_20, args.out_dir)

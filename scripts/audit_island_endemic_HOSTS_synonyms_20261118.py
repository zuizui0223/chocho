#!/usr/bin/env python3
"""Independent island butterfly HOSTS coverage after guarded GBIF taxon reconciliation.

The biological denominator is the source-frozen INBO island checklist, not HOSTS.
Accepted name/synonym matches are DATABASE COVERAGE only, not field host absence
or verified establishment. Taxonomic disagreements and GBIF outages remain UNKNOWN.
"""
import argparse
import csv
import datetime as dt
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

BASE = Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
FROZEN = BASE / "independent_42_species_vs_frozen239.csv"
HOSTS_SHA = "808e0b869f9ec1adf8efff87cf6a395adda103e0"
ORIGINAL_CSV_BLOB = "3c56054fa21e667eed422a67010059c8e198b1bc"
API = "https://api.gbif.org/v1/species"

# Inclusion is independently justified, not selected by HOSTS records.
# Canary Island endemic taxa are one related evolutionary radiation; they
# are NOT five independent island-origin/phylogeny replications.
COHORTS = {
    "Madeira_native_nonendemic": {
        "source": "Borges et al. 2008 pp342-343",
        "origin_class": "native_nonendemic",
        "names": [
            "Colias crocea", "Danaus plexippus", "Lampides boeticus",
            "Leptotes pirithous", "Pararge aegeria", "Vanessa atalanta",
            "Vanessa cardui",
        ],
    },
    "Madeira_extant_endemic": {
        "source": "Borges et al. 2008 pp342-343",
        "origin_class": "island_endemic_species",
        "names": [
            "Gonepteryx maderensis", "Hipparchia maderensis", "Pararge xiphia"
        ],
    },
    "Canary_single_island_endemic_Hipparchia": {
        "source": "Spilani et al. 2025, Insect Systematics and Diversity, doi:10.1093/isd/ixaf043",
        "origin_class": "island_endemic_species",
        "names": [
            "Hipparchia tamadabae", "Hipparchia wyssii", "Hipparchia gomera",
            "Hipparchia tilosi", "Hipparchia bacchus"
        ],
    },
    "Azores_endemic_Hipparchia": {
        "source": "Azores regional authority SIARAM arthropod conservation inventory (Hipparchia azorina, H. miguelensis)",
        "origin_class": "island_endemic_species",
        "names": ["Hipparchia azorina", "Hipparchia miguelensis"],
    },
}


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def canonical_name(raw):
    parts = (raw or "").strip().split()
    if len(parts) != 2:
        return ""
    if not re.fullmatch(r"[A-Za-z][A-Za-z-]*", parts[0]):
        return ""
    if not re.fullmatch(r"[a-z][a-z-]*", parts[1]):
        return ""
    return (parts[0] + " " + parts[1]).lower()


def host_has_genus(raw):
    genus = (raw or "").strip()
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z-]*", genus)) and genus.lower() not in {
        "unknown", "unidentified", "indet", "fungi"
    }


def request_json(url, receipt_path, offline=False):
    if receipt_path.exists():
        return json.loads(receipt_path.read_text(encoding="utf-8"))
    if offline:
        raise RuntimeError("Missing source receipt in offline mode: " + str(receipt_path))
    last = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "chocho-source-locked-independent-host-coverage/2026-11-18 (research audit)"
            })
            with urllib.request.urlopen(request, timeout=25) as response:
                body = json.load(response)
            receipt_path.parent.mkdir(parents=True, exist_ok=True)
            receipt_path.write_text(json.dumps({
                "url": url,
                "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "response": body
            }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            return {"response": body}
        except (OSError, ValueError, urllib.error.URLError) as exc:
            last = str(exc)
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(last)


def parse_api(receipt):
    return receipt.get("response", {})


def obtain(url, path, offline=False):
    return parse_api(request_json(url, path, offline=offline))


def make_hosts_index(hosts_csv):
    index = defaultdict(lambda: {"rows": 0, "named_genus": 0, "named_species": 0})
    source_n = 0
    with hosts_csv.open(encoding="utf-8-sig", newline="") as handle:
        for item in csv.DictReader(handle):
            source_n += 1
            name = canonical_name(
                ((item.get("Insect Genus") or "") + " " +
                 (item.get("Insect Species") or "")).strip())
            family = (item.get("Insect Family") or "").strip().lower()
            if not name or not family:
                continue
            record = index[(name, family)]
            record["rows"] += 1
            if host_has_genus(item.get("Hostplant Genus", "")):
                record["named_genus"] += 1
                epithet = (item.get("Hostplant Species") or "").strip()
                if re.fullmatch(r"[a-z][a-z-]*", epithet) and epithet.lower() not in {
                    "sp", "spp", "indet", "unknown", "cf", "aff"
                }:
                    record["named_species"] += 1
    if source_n != 140485:
        raise ValueError("Pinned HOSTS rows unexpected; stop: " + str(source_n))
    return index


def gbif_name_set(name, family, receipts, offline):
    safe = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    query = urllib.parse.urlencode({
        "name": name, "family": family, "order": "Lepidoptera",
        "kingdom": "Animalia", "strict": "true",
    })
    match_url = API + "/match?" + query
    result = obtain(match_url, receipts / (safe + "_match.json"), offline)
    guard = {
        "match_type": result.get("matchType", "NONE"),
        "confidence": result.get("confidence"),
        "rank": result.get("rank"),
        "matched_name": result.get("scientificName", ""),
        "matched_family": result.get("family", ""),
        "status": result.get("status", "")
    }
    if (guard["match_type"] != "EXACT" or
        (guard["confidence"] or 0) < 90 or
        guard["rank"] != "SPECIES" or
        guard["matched_family"].lower() != family.lower() or
        guard["status"] == "DOUBTFUL"):
        return None, set(), guard, "gbif_match_guard_failed"
    accepted_key = (result.get("acceptedUsageKey") or result.get("acceptedKey") or
                    result.get("speciesKey") or result.get("usageKey"))
    if not accepted_key:
        return None, set(), guard, "no_accepted_key"
    detail = obtain(API + "/" + str(accepted_key),
                    receipts / (safe + "_accepted.json"), offline)
    if (detail.get("rank") != "SPECIES" or
        (detail.get("family") or "").lower() != family.lower() or
        detail.get("taxonomicStatus", "ACCEPTED") not in ("ACCEPTED", "PROVISIONALLY_ACCEPTED")):
        return None, set(), guard, "accepted_species_guard_failed"
    accepted_name = canonical_name(detail.get("canonicalName", ""))
    if not accepted_name:
        return None, set(), guard, "no_accepted_species_binomial"
    synonyms = []
    offset = 0
    while True:
        url = API + "/" + str(accepted_key) + "/synonyms?" + urllib.parse.urlencode({
            "limit": 500, "offset": offset
        })
        page = obtain(url, receipts / (safe + "_synonyms_" + str(offset) + ".json"), offline)
        synonyms.extend(page.get("results", []))
        if page.get("endOfRecords", True) or len(page.get("results", [])) == 0:
            break
        offset += len(page["results"])
        if offset > 4000:
            return None, set(), guard, "too_many_synonyms_for_safe_audit"
    candidates = {canonical_name(name), accepted_name}
    for synonym in synonyms:
        if synonym.get("rank") != "SPECIES":
            continue
        syn_family = synonym.get("family")
        if syn_family and syn_family.lower() != family.lower():
            continue
        syn = canonical_name(synonym.get("canonicalName", ""))
        if syn:
            candidates.add(syn)
    return str(accepted_key), candidates, guard, "resolved"


def analyze(frozen_csv, hosts_csv, output, offline=False):
    output.mkdir(parents=True, exist_ok=True)
    receipts = output / "gbif_source_receipts"
    frozen = rows(frozen_csv)
    if len(frozen) != 42 or sum(int(r["hosts_named_genus"]) for r in frozen) != 19:
        raise ValueError("Frozen independent denominator changed")
    if sum(int(r["in_frozen239"]) for r in frozen) != 12:
        raise ValueError("Frozen original GEB selected cohort changed")
    index = make_hosts_index(hosts_csv)
    by_name = {r["scientific_name"]: r for r in frozen}
    for cohort_name, cohort in COHORTS.items():
        if len(set(cohort["names"])) != len(cohort["names"]):
            raise ValueError("Duplicate source taxon in " + cohort_name)
        if any(x not in by_name for x in cohort["names"]):
            raise ValueError("Source-cohort taxon missing from frozen INBO 42: " + cohort_name)
    records = []
    for i, source in enumerate(frozen):
        name = source["scientific_name"]
        family = source["checklist_family"]
        original = canonical_name(name)
        exact = int(source["hosts_named_genus"])
        try:
            key, names, guard, status = gbif_name_set(name, family, receipts, offline)
        except Exception as exc:
            key, names, guard, status = None, set(), {}, "gbif_api_unavailable"
            guard["error"] = str(exc)
        new = sorted(n for n in names if n != original and index.get(
            (n, family.lower()), {}).get("named_genus", 0) > 0)
        row = {
            "scientific_name": name, "family": family,
            "hosts_exact_genus": exact,
            "hosts_exact_species": int(source["hosts_named_species"]),
            "in_original239": int(source["in_frozen239"]),
            "gbif_status": status,
            "gbif_accepted_species_key": key or "",
            "gbif_match_guard": guard,
            "gbif_synonym_n": max(0, len(names) - 1),
            "candidate_alias_host_binomials": new,
            "hosts_genus_after_guarded_synonyms": int(bool(exact or new)) if status == "resolved" else exact,
            "additional_match": bool(not exact and new and status == "resolved"),
            "unresolved_potential_additional": bool(not exact and status != "resolved"),
        }
        records.append(row)
    # Crucial: prevent accepted species collisions from silently collapsing
    # distinct checklist taxa (possible in island species complexes).
    key_groups = defaultdict(list)
    for record in records:
        if record["gbif_accepted_species_key"]:
            key_groups[record["gbif_accepted_species_key"]].append(record)
    collisions = []
    for key, members in key_groups.items():
        if len(members) <= 1:
            continue
        collisions.append({"accepted_key": key,
                           "distinct_checklist_taxa": [r["scientific_name"] for r in members]})
        for member in members:
            member["gbif_status"] = "ambiguous_accepted_species_collision"
            member["hosts_genus_after_guarded_synonyms"] = member["hosts_exact_genus"]
            member["additional_match"] = False
            member["unresolved_potential_additional"] = not bool(member["hosts_exact_genus"])
    by_name = {r["scientific_name"]: r for r in records}
    cohort_results = {}
    for cohort_name, cohort in COHORTS.items():
        subset = [by_name[n] for n in cohort["names"]]
        cohort_results[cohort_name] = {
            "source": cohort["source"],
            "origin_class": cohort["origin_class"],
            "n": len(subset),
            "exact_hosts_genus": sum(r["hosts_exact_genus"] for r in subset),
            "confirmed_name_reconciled_hosts_genus": sum(r["hosts_genus_after_guarded_synonyms"] for r in subset),
            "gbif_unresolved": sum(r["gbif_status"] != "resolved" for r in subset),
            "in_original239": sum(r["in_original239"] for r in subset),
            "taxa": cohort["names"],
        }
    extras = [r["scientific_name"] for r in records if r["additional_match"]]
    pending = [r["scientific_name"] for r in records if r["unresolved_potential_additional"]]
    summary = {
        "schema": "chocho_independent_island_endemics_guarded_taxonomy_v1",
        "analysis_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "taxonomic_source": "GBIF v1 backbone API, response snapshots stored",
        "host_source_commit": HOSTS_SHA,
        "original_42_csv_git_blob": ORIGINAL_CSV_BLOB,
        "meaning": "reported HOSTS larval plant genus coverage after guarded accepted-name/synonym matching; not biological host completeness",
        "independent_unique_butterflies": len(records),
        "original_exact_binomial_named_host": sum(r["hosts_exact_genus"] for r in records),
        "guarded_synonym_named_host_confirmed_min": sum(r["hosts_genus_after_guarded_synonyms"] for r in records),
        "guarded_additions": extras,
        "unresolved_no_exact_host": pending,
        "all_42_gbif_fully_resolved": all(r["gbif_status"] == "resolved" for r in records),
        "gbif_taxon_collisions": collisions,
        "cohorts": cohort_results,
        "limitations": [
            "Alias reconciliation is a taxonomic matching audit of HOSTS, not a biological absence estimate.",
            "Fuzzy, rank-inconsistent and multi-species GBIF matches are not automatically resolved.",
            "Accepted-name collision quarantined; island endemics not silently merged into mainland species.",
            "Madeira vs Canary vs Azores are differing independent source scopes and dates.",
            "Five Canary Hipparchia are one radiation, not independent species-richness replicates.",
            "Even a full taxonomy match cannot retrieve hosts never entered in NHM HOSTS."
        ]
    }
    with (output / "accepted_name_host_coverage_42.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        cols = [x for x in records[0] if x != "gbif_match_guard"]
        writer = csv.DictWriter(handle, fieldnames=cols)
        writer.writeheader()
        for record in records:
            row = {x: record[x] for x in cols}
            row["candidate_alias_host_binomials"] = "; ".join(row["candidate_alias_host_binomials"])
            writer.writerow(row)
    (output / "accepted_name_host_coverage_42.json").write_text(
        json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (output / "synonym_origin_generalization_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("HOSTS_ALIAS_AUDIT_SUMMARY=" + json.dumps({
        "exact_42": summary["original_exact_binomial_named_host"],
        "guarded_reconciled_min_42": summary["guarded_synonym_named_host_confirmed_min"],
        "added": extras,
        "unresolved_no_exact": pending,
        "gbif_fully_resolved": summary["all_42_gbif_fully_resolved"],
        "cohorts": {name: {
            "n": v["n"], "exact": v["exact_hosts_genus"],
            "reconciled_min": v["confirmed_name_reconciled_hosts_genus"],
            "unresolved": v["gbif_unresolved"], "in_239": v["in_original239"]
        } for name, v in cohort_results.items()}
    }, ensure_ascii=False))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-42", type=Path, default=FROZEN)
    parser.add_argument("--hosts", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    analyze(args.frozen_42, args.hosts, args.out_dir, args.offline)

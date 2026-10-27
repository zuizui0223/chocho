#!/usr/bin/env python3
"""Original HOSTS Location x WCVP accepted-host plant native/introduced audit.

Source-only ecological record geography, NOT a measured butterfly host
switch, dispersal, establishment, colonization, fitness or causal impact.
Deduplicate to butterfly x accepted host x original Location; preserve
original raw records only as aggregate counts and NEVER export host edges.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from build_hosts_location_scope_crosswalk_20261027 import classify

SCHEMA = "chocho_original_hosts_location_wcvp_plant_origin_20261027_v01"
SOURCE_ZIP_SHA = "408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
COUNTRY_CROSSWALK_COMMIT = "4104c784ddcf97318c4a1bb8bd949b3147628494"
EXPECTED_ATOM_SHA = "050294627df4b0dd0dbfe635888cd2396ddbdd75acfc69b811f9f03c03a35370"
EXPECTED_HOSTS_ROWS = 140485
EXPECTED_SPECIES_HOST_ROWS = 86826
EXPECTED_LOCATION_LABELS = 257
# These are direct, independently-checkable botanical region names,
# not assumed classifications of their broader parent countries.
SINGLE_L3_ALIASES = {
    "Hawaii": ("US", "HAW"),
    "West Malaysia": ("MY", "MLY"),
    # East Malaysia, Sabah, Sarawak, Java, Ryukyu are NOT reduced to one
    # botanical unit without a separate verified source geometry.
}
LABELS_DIAGNOSTIC = (
    "USA", "Brazil", "Taiwan", "Japan", "India", "West Malaysia",
    "Hawaii", "Nearctic", "Oriental", "USA [interception]",
)
STATUS = (
    "NATIVE_CONFIRMED", "INTRODUCED_CONFIRMED", "WCVP_UNLISTED",
    "MIXED_OR_PARTIAL", "COUNTRY_MAPPING_UNSAFE",
    "SUBREGION_MAPPING_UNSAFE", "UNMAPPED_OR_UNREVIEWED_LOCATION",
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def accepted_pairs(path, focal):
    ret = defaultdict(set)
    families = defaultdict(set)
    for row in read_rows(path):
        if row["insect_species"] not in focal:
            continue
        species = row["insect_species"]
        raw = row["input_host_name"]
        id_ = str(row["accepted_plant_name_id"]).strip()
        if id_:
            ret[species, raw].add(id_)
            families[id_].add(row.get("family", "").strip() or "<UNKNOWN>")
    ambiguous = sum(len(v) > 1 for v in ret.values())
    return {k: next(iter(ids)) for k, ids in ret.items() if len(ids) == 1}, {
        "source_distinct_pair_names": len(ret),
        "source_ambiguous_accepted_id_names": ambiguous,
        "source_distinct_plant_ids": len(set().union(*ret.values())) if ret else 0,
    }, {k: sorted(v) for k, v in families.items()}


def wcvp_range(native_path, current_path):
    def load(path):
        table = defaultdict(set)
        for row in read_rows(path):
            pid = str(row["accepted_plant_name_id"]).strip()
            reg = str(row["area_code_l3"]).strip()
            if pid and reg:
                table[pid].add(reg)
        return table
    native, current = load(native_path), load(current_path)
    discrepancy = {p: len(locs-current.get(p, set())) for p, locs in native.items()
                   if not locs <= current.get(p, set())}
    if discrepancy:
        raise ValueError("Original WCVP native ranges are not subset of current "
                         + str(len(discrepancy)))
    return native, current


def load_country_atoms(atoms_path, atoms_manifest_path, countries_path):
    manifest = json.loads(Path(atoms_manifest_path).read_text())
    summary = json.loads(Path(countries_path).read_text())
    if manifest["atoms_csv_sha256"] != sha256(atoms_path) or \
       manifest["atoms_csv_sha256"] != EXPECTED_ATOM_SHA:
        raise ValueError("Country source disjoint atom SHA mismatch")
    if manifest.get("crosswalk_summary_sha256") != sha256(countries_path):
        raise ValueError("Original country crosswalk SHA mismatch")
    if summary.get("schema") != "country_wgsrpd_area_crosswalk_v1":
        raise ValueError("Not original WCVP/TDWG country crosswalk")
    if manifest["atom_rows"] != 969 or summary["country_count"] != 247:
        raise ValueError("Area atoms or political country roster has changed")
    atom_meta = {r["country"]: r for r in manifest["countries"]}
    by_country = defaultdict(list)
    for row in read_rows(atoms_path):
        country, parts = row["country"], row["wgsrpd3_members"].split("|")
        area = float(row["area_km2"])
        if not parts or len(set(parts)) != len(parts) or area <= 0:
            raise ValueError("Invalid country disjoint atom membership")
        by_country[country].append((tuple(parts), area))
    if set(by_country) != set(atom_meta):
        raise ValueError("Some countries have zero or unverified geographic atoms")
    area_total = 0
    regions = set()
    for code, entries in by_country.items():
        m = atom_meta[code]
        total = sum(a for _, a in entries)
        if not math.isclose(total, float(m["atoms_area_km2"]),
                            rel_tol=1e-8, abs_tol=0.01):
            raise ValueError("Country atom subarea does not reproduce original source")
        if not math.isclose(total/float(m["country_area_km2"]),
                            float(m["guard_coverage"]), rel_tol=1e-6, abs_tol=1e-6):
            raise ValueError("Country atom coverage drift")
        if m["area_identity_checks"] != "PASS":
            raise ValueError("Original country geometry identity check failed")
        area_total += total
        for parts, _ in entries:
            regions.update(parts)
    return by_country, atom_meta, regions, {
        "source_country_crosswalk_commit": COUNTRY_CROSSWALK_COMMIT,
        "area_atoms_sha256": sha256(atoms_path),
        "area_atoms_summary_sha256": sha256(atoms_manifest_path),
        "crosswalk_summary_sha256": sha256(countries_path),
        "country_count": len(by_country),
        "atom_count": sum(map(len, by_country.values())),
        "distinct_original_wgsrpd3_regions": len(regions),
        "total_mapped_area_km2": area_total,
    }


def per_region(pid, region, native, current):
    if region in native.get(pid, ()):
        return "NATIVE_CONFIRMED"
    if region in current.get(pid, ()):
        return "INTRODUCED_CONFIRMED"
    return "WCVP_UNLISTED"


def dominant_country(pid, country, *, threshold, atoms, metadata, native, current):
    """Area-dominance over WHOLE country, not botanical-region counts."""
    m = metadata.get(country)
    if not m or country not in atoms or \
       float(m["guard_coverage"]) < .95-1e-9 or \
       not m.get("area_identity_checks") == "PASS":
        return "COUNTRY_MAPPING_UNSAFE"
    area_country = float(m["country_area_km2"])
    if not area_country > 0:
        return "COUNTRY_MAPPING_UNSAFE"
    amounts = Counter()
    for members, area in atoms[country]:
        types = {per_region(pid, r, native, current) for r in members}
        if len(types) == 1:
            amounts[next(iter(types))] += area
        else:
            amounts["ambiguous_botanical_overlap"] += area
    # Country uncovered area is NOT a WCVP-UNLISTED plant.
    amounts["unmapped_country_area"] = max(0., area_country - sum(
        area for members, area in atoms[country]))
    for label in ("NATIVE_CONFIRMED", "INTRODUCED_CONFIRMED", "WCVP_UNLISTED"):
        if amounts[label] / area_country >= threshold - 1e-10:
            return label
    return "MIXED_OR_PARTIAL"


def subregion_status(pid, country, region, *, atoms, original_regions,
                     native, current):
    if region not in original_regions or country not in atoms or \
       not any(region in parts for parts, _ in atoms[country]):
        return "SUBREGION_MAPPING_UNSAFE"
    return per_region(pid, region, native, current)


def classify_location(pid, label, *, threshold, atoms, atom_metadata,
                      original_regions, native, current):
    scope, country = classify(label)
    if scope == "SINGLE_COUNTRY_OR_TERRITORY_NAME":
        return dominant_country(pid, country, threshold=threshold,
                                atoms=atoms, metadata=atom_metadata,
                                native=native, current=current)
    if scope == "SUBNATIONAL_NAME_NEEDS_GEOMETRY":
        if label not in SINGLE_L3_ALIASES:
            return "SUBREGION_MAPPING_UNSAFE"
        nation, region = SINGLE_L3_ALIASES[label]
        if nation != country:
            return "SUBREGION_MAPPING_UNSAFE"
        return subregion_status(pid, nation, region, atoms=atoms,
                                original_regions=original_regions,
                                native=native, current=current)
    return "UNMAPPED_OR_UNREVIEWED_LOCATION"


def observations(hosts_path, accepted_map, focal, *, native, current, atoms,
                 atom_metadata, source_regions, families):
    counters = Counter()
    triples = set()
    raw_counts = Counter()
    source_labels = Counter()
    focal_locations = defaultdict(Counter)
    for row in read_rows(hosts_path):
        counters["source_all_raw_rows"] += 1
        fields=("Insect Genus", "Insect Species", "Hostplant Genus", "Hostplant Species")
        if not all(str(row.get(c) or "").strip() for c in fields):
            continue
        counters["source_full_species_binomial_rows"] += 1
        species = (row["Insect Genus"].strip() + " " +
                   row["Insect Species"].strip())
        if species not in focal:
            continue
        counters["source_focal_butterfly_raw_species_rows"] += 1
        host = row["Hostplant Genus"].strip()+" "+row["Hostplant Species"].strip()
        pid = accepted_map.get((species, host))
        if pid is None:
            counters["source_focal_rows_unmatched_exact_accepted_pair"] += 1
            continue
        counters["source_focal_accepted_pair_raw_rows"] += 1
        label = str(row.get("Location") or "").strip()
        if not label:
            counters["source_focal_accepted_rows_empty_location"] += 1
        raw_counts[species] += 1
        source_labels[label] += 1
        focal_locations[species][label] += 1
        triples.add((species, pid, label))
    if counters["source_all_raw_rows"] != EXPECTED_HOSTS_ROWS or \
       counters["source_full_species_binomial_rows"] != EXPECTED_SPECIES_HOST_ROWS:
        raise ValueError("Pinned raw HOSTS source row accounting changed")
    if not triples:
        raise ValueError("Original accepted HOSTS locations have no 239-butterfly support")

    status_lookup = {}
    triple_status = defaultdict(list)
    report = {}
    for cutoff in (.95, .8):
        by_status = Counter()
        by_raw = Counter()
        links = defaultdict(set)
        pairs_native = set()
        all_triples = {}
        source_label_evidence = defaultdict(Counter)
        for sp, pid, label in sorted(triples):
            cache_key = (pid, label, cutoff)
            if cache_key not in status_lookup:
                status_lookup[cache_key] = classify_location(
                    pid, label, threshold=cutoff, atoms=atoms,
                    atom_metadata=atom_metadata, original_regions=source_regions,
                    native=native, current=current)
                if status_lookup[cache_key] not in STATUS:
                    raise ValueError("Unexpected plant origin status")
            status = status_lookup[cache_key]
            all_triples[sp, pid, label] = status
            if status == "NATIVE_CONFIRMED":
                pairs_native.add((sp, pid))
            by_status[status] += 1
            links[status].add((sp, pid))
            source_label_evidence[label][status] += 1
        intro = [pair for pair, status in all_triples.items()
                 if status == "INTRODUCED_CONFIRMED"]
        intro_pairs = {(sp, pid) for sp, pid, loc in intro}
        has_native_intro = [((sp,pid) in pairs_native) for sp,pid,loc in intro]
        matched_intro_pairs = intro_pairs & pairs_native
        sampled_native = Counter(sp for sp,pid in pairs_native)
        sampled_intro = Counter(sp for sp,pid in intro_pairs)
        for r in read_rows(hosts_path):
            sp = ((r.get("Insect Genus") or "").strip() + " " +
                  (r.get("Insect Species") or "").strip())
            if sp not in focal:
                continue
            host = ((r.get("Hostplant Genus") or "").strip() + " " +
                    (r.get("Hostplant Species") or "").strip())
            pid = accepted_map.get((sp,host))
            if pid is None:
                continue
            label = str(r.get("Location") or "").strip()
            status = all_triples[sp,pid,label]
            by_raw[status]+=1
        by_country = {}
        for loc in LABELS_DIAGNOSTIC:
            source_scope, _ = classify(loc)
            by_country[loc] = {
                "original_report_rows_with_239_matched_host_pair": focal_locations_count(
                    loc, focal_locations
                ),
                "geographic_scope":source_scope,
                "unique_triples_by_plant_origin":dict(source_label_evidence.get(loc,{})),
            }
        family_intro = Counter()
        for sp,pid in intro_pairs:
            possible=families.get(pid, [])
            family_intro[possible[0] if len(possible)==1 else "<FAMILY_UNCERTAIN>"]+=1
        # Statistics are aggregate, NOT organismal sample size or source-independent events.
        report[str(cutoff)] = {
            "unique_species_host_location_triples":len(triples),
            "triples_by_plant_origin":dict(by_status),
            "original_raw_records_by_plant_origin":dict(by_raw),
            "distinct_species_host_pairs_by_plant_origin":{
                s:len(v) for s,v in links.items()},
            "confirmed_introduced_plant_area_triples":len(intro),
            "introduced_triples_with_same_pair_documented_in_native_area":sum(
                has_native_intro),
            "introduced_triples_without_native_area_pair_report":sum(
                not v for v in has_native_intro),
            "introduced_distinct_species_host_pairs":len(intro_pairs),
            "introduced_pairs_also_documented_in_native_area":len(matched_intro_pairs),
            "introduced_pairs_without_native_area_report":len(
                intro_pairs - pairs_native),
            "different_butterfly_species_with_introduced_plant_area_report":len(
                sampled_intro),
            "different_butterfly_species_with_native_plant_area_report":len(
                sampled_native),
            "introduced_host_report_species_with_any_native_pair_evidence":sum(
                sp in sampled_native for sp in sampled_intro),
            "top_15_plant_families_for_introduced_source_pairs": [
                {"family": name, "pair_count": n}
                for name,n in family_intro.most_common(15)
            ],
            "predeclared_source_labels":by_country,
            "observational_warning": (
                "Introduced-record-only host pairs are NOT confirmed newly "
                "adopted larval hosts. Absent native HOSTS citation is missing "
                "documentation, not ecological native host non-use."
            ),
        }
    return {
        "record_flow":dict(counters),
        "unique_focal_species_host_location_triples":len(triples),
        "source_distinct_location_labels_for_focal_taxa":len(source_labels),
        "distinct_focal_butterfly_species_in_source":len(raw_counts),
        "primary_and_sensitivity":report,
    }


def focal_locations_count(loc, by_species):
    return sum(s[loc] for s in by_species.values())


def run(args):
    original = json.loads((args.author_source / "SOURCE_MANIFEST.json").read_text())
    if original.get("inputs_count") != 11:
        raise ValueError("Author source manifest drifted")
    for path, expected in original["sha256"].items():
        if sha256(args.author_source / path) != expected:
            raise ValueError("Source WCVP/HOSTS underlying member hash mismatch "+path)
    root=args.author_source/"inputs"
    roster=read_rows(root/"anthropogenic_species_metrics.csv")
    focal={r["species"] for r in roster}
    if len(focal)!=239 or len(roster)!=239:
        raise ValueError("Not the original 239 butterfly source cohort")
    accepted, accepted_meta, families=accepted_pairs(
        root/"original/insect_host_accepted.csv",focal)
    if accepted_meta["source_ambiguous_accepted_id_names"]:
        raise ValueError("Original source accepted host pair names ambiguous")
    native,current=wcvp_range(root/"original/native_extant_nondoubtful_wgsrpd3.csv",
                              root/"original/extant_nondoubtful_wgsrpd3.csv")
    atoms, atom_meta, regions, geography=load_country_atoms(
        args.country_atoms,args.country_atoms_manifest,args.country_summary)
    rows=observations(args.hosts_resource,accepted,focal,native=native,
                      current=current,atoms=atoms,atom_metadata=atom_meta,
                      source_regions=regions,families=families)
    obj={
        "schema":SCHEMA,
        "status":"SOURCE_LOCKED_RETROSPECTIVE_ORIGIN_RECORDS_NO_HOST_ADOPTION_CLAIM",
        "source": {
            "HOSTS_mirror":"globalbioticinteractions/HOSTS@808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "HOSTS_resource_sha256":sha256(args.hosts_resource),
            "author_only_inner_archive_sha256":SOURCE_ZIP_SHA,
            "accepted_pairs_sha256":sha256(root/"original/insect_host_accepted.csv"),
            "native_WCVP_sha256":sha256(root/"original/native_extant_nondoubtful_wgsrpd3.csv"),
            "current_WCVP_sha256":sha256(root/"original/extant_nondoubtful_wgsrpd3.csv"),
            "original_239_source_sha256":sha256(root/"anthropogenic_species_metrics.csv"),
            "original_accepted_join":accepted_meta,
            "source_native_accepted_plant_ids":len(native),
            "source_current_accepted_plant_ids":len(current),
            "geography":geography,
        },
        "result":rows,
        "inference_limits":[
            "This is a report-location botanical origin classification, not a measured insect invasion or larval host-switch event.",
            "No native-area record of a butterfly-host combination is not proof it did not occur there.",
            "Introduced accepted WCVP plant distribution can include heterogeneous source landscapes and dates, and does not prove same-site host plant availability.",
            "WCVP unlisted region is not cultivated-only, confirmed absence, or ecological treatment control.",
            "Source butterfly native/introduced geography not established by the plant's botanical status.",
            "The 239 butterflies are not a population-independent sample; publication and source recording effort are highly uneven.",
            "Country-level source area-weighted status and exact verified subnational TDWG3 status are distinct spatial grains.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    main=obj["result"]["primary_and_sensitivity"]["0.95"]
    print(json.dumps({
        "source_239_butterflies":239,
        "source_counts":obj["result"]["record_flow"],
        "source_distinct_triples":main["unique_species_host_location_triples"],
        "main95":{k:main[k] for k in (
            "triples_by_plant_origin","confirmed_introduced_plant_area_triples",
            "introduced_triples_with_same_pair_documented_in_native_area",
            "introduced_triples_without_native_area_pair_report",
            "introduced_distinct_species_host_pairs",
            "introduced_pairs_also_documented_in_native_area",
            "different_butterfly_species_with_introduced_plant_area_report")},
        "sensitivity80":obj["result"]["primary_and_sensitivity"]["0.8"][
            "triples_by_plant_origin"
        ],
    },ensure_ascii=False,sort_keys=True))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",type=Path,required=True)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--country-atoms",type=Path,required=True)
    p.add_argument("--country-atoms-manifest",type=Path,required=True)
    p.add_argument("--country-summary",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())

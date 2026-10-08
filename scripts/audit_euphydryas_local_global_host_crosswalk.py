#!/usr/bin/env python3
"""Input-data audit, NOT an ecological survival or removal experiment.

Reconciles named local E. editha hosts with frozen species-wide HOSTS x WCVP
relations and tests only static uniqueness of associated introduced regions.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def csv_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f)


def source_hosts(path, consumer, focal):
    matches = {p: 0 for p in focal}
    seen_consumer = 0
    for r in csv_rows(path):
        insect = " ".join([
            (r.get("Insect Genus") or "").strip(),
            (r.get("Insect Species") or "").strip()
        ]).strip()
        if insect != consumer:
            continue
        seen_consumer += 1
        host = " ".join([
            (r.get("Hostplant Genus") or "").strip(),
            (r.get("Hostplant Species") or "").strip()
        ]).strip()
        if host in matches:
            matches[host] += 1
    return matches, seen_consumer


def distribution(path):
    out = defaultdict(set)
    for r in csv_rows(path):
        host = (r.get("accepted_plant_name_id") or "").strip()
        region = (r.get("area_code_l3") or "").strip()
        if host and region:
            out[host].add(region)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol-json", required=True, type=Path)
    parser.add_argument("--original-hosts-csv", required=True, type=Path)
    parser.add_argument("--interaction-csv", required=True, type=Path)
    parser.add_argument("--native-csv", required=True, type=Path)
    parser.add_argument("--contemporary-csv", required=True, type=Path)
    parser.add_argument("--global-descriptors-csv", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    spec = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if spec["schema"] != "chocho_euphydryas_four_hosts_exact_crosswalk_v0.1":
        raise RuntimeError("Unknown frozen crosswalk protocol schema")
    butterfly = spec["species"]
    focal = [p["binomial"] for p in spec["focal_plants"]]
    if len(focal) != len(set(focal)) or len(focal) != 4:
        raise RuntimeError("Four distinct a priori local plants expected")

    raw_match, raw_consumer_rows = source_hosts(args.original_hosts_csv, butterfly, focal)
    # Use accepted taxon IDs, not host genus approximations.
    all_ids = set()
    exact_links = defaultdict(list)
    all_links = []
    for row in csv_rows(args.interaction_csv):
        if (row.get("insect_species") or "").strip() != butterfly:
            continue
        hid = (row.get("accepted_plant_name_id") or "").strip()
        if not hid:
            raise RuntimeError("Missing accepted taxon ID")
        all_ids.add(hid)
        source_name = (row.get("input_host_name") or "").strip()
        accepted_name = (row.get("accepted_name") or "").strip()
        all_links.append((source_name, accepted_name, hid))
        for p in focal:
            if p in (source_name, accepted_name):
                exact_links[p].append({
                    "input_host_name": source_name,
                    "accepted_name": accepted_name,
                    "accepted_plant_name_id": hid
                })
    native = distribution(args.native_csv)
    contemporary = distribution(args.contemporary_csv)
    for hid in all_ids:
        if native[hid] - contemporary[hid]:
            raise RuntimeError(f"Native region not present contemporary for accepted ID {hid}")

    N = set().union(*(native[hid] for hid in all_ids))
    C = set().union(*(contemporary[hid] for hid in all_ids))
    G = C - N

    descriptors = [r for r in csv_rows(args.global_descriptors_csv)
                   if (r.get("species") or "").strip() == butterfly]
    if len(descriptors) != 1:
        raise RuntimeError("Global descriptor must match exact butterfly once")
    d = descriptors[0]
    actual = {
        "resolved_host_species": len(all_ids),
        "native_resource_units": len(N),
        "contemporary_resource_units": len(C),
        "introduced_added_units": len(G)
    }
    expected = spec["fixed_expected_global"]
    if actual != expected or any(int(d[k]) != actual[k] for k in actual):
        raise RuntimeError(
            "Crosswalk input is not the frozen published E. editha input: "
            + json.dumps({"expected": expected, "actual": actual,
                          "descriptor": {k: d[k] for k in actual}})
        )

    observations = []
    for plant in spec["focal_plants"]:
        p = plant["binomial"]
        links = exact_links[p]
        ids = sorted({item["accepted_plant_name_id"] for item in links})
        rec = {
            "plant": p,
            "published_local_role": plant["role"],
            "source_publications": plant["references"],
            "raw_Hosts_exact_butterfly_link_rows": raw_match[p],
            "WCVP_resolved_exact_accepted_taxon_ids": ids,
            "accepted_link_details": links,
            "link_state": (
                "NO_EXACT_HOSTS_LINK" if raw_match[p] == 0 and not ids else
                "RAW_LINK_UNRESOLVED_BY_WCVP" if raw_match[p] and not ids else
                "EXACT_ACCEPTED_ID_LINK" if len(ids) == 1 else
                "ACCEPTED_TAXON_AMBIGUITY_OR_MULTIPLE_IDENTITIES"
            )
        }
        if len(ids) == 1:
            hid = ids[0]
            mapped_n = native[hid]
            mapped_c = contemporary[hid]
            introduced = mapped_c - mapped_n
            # Counterfactual: retain focal host's native distribution, leave every
            # other *fixed recorded* host untouched. This is not botanical removal.
            after = set(mapped_n)
            after.update(*(contemporary[other] for other in all_ids if other != hid))
            lost = C - after
            rec.update({
                "botanical_native_region_count": len(mapped_n),
                "botanical_contemporary_region_count": len(mapped_c),
                "botanical_introduced_only_region_count": len(introduced),
                "focal_butterfly_sole_added_resource_regions": len(lost),
                "focal_butterfly_sole_added_resource_fraction": (
                    len(lost) / len(G) if G else None
                ),
                "focal_introduced_botanical_regions_already_covered_by_other_hosts":
                    len(introduced - lost),
                "missing_only_by_removing_focal_introduced_distribution":
                    len(lost),
            })
            if not lost.issubset(G):
                raise RuntimeError("Host introduced-only removal must not destroy native opportunity")
        observations.append(rec)

    result = {
        "schema": "chocho_euphydryas_four_local_hosts_exact_global_crosswalk_result_v0.1",
        "status": "POSTHOC_EXACT_INPUT_AUDIT_NOT_LOCAL_FITNESS",
        "input_provenance": spec["source"],
        "species": butterfly,
        "global_reconstruction_checks": actual,
        "original_HOSTS_insect_rows": raw_consumer_rows,
        "local_literature_hosts": observations,
        "limits": [
            "Missing recorded link does not mean a butterfly cannot use the plant",
            "Species-wide host relations and coarse botanical ranges do not measure meadow-specific host use",
            "Region counts are NOT local population fitness or demographic fallback",
            "Static host deletion is not a plant-removal experiment",
            "Four examples selected from published cases are post-hoc and not a random sample",
        ]
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

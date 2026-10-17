#!/usr/bin/env python3
"""Paired geographic-identity turnover under two larval-host inventories.

The comparison is between the SAME butterfly x WGSRPD3 candidates under the
archived 239-butterfly HOSTS and 81-species BCE/Clarke host-augmented input.
It is a descriptive DATA-SOURCE sensitivity, not a temporal transition or
statistical host-choice mechanism. No native/introduced geographic identities
are fabricated; pairwise intersection/union is algebraically recovered from
the previously audited transition counts.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

SCHEMA = "chocho_butterfly_bce_host_inventory_resource_novelty_v0.1"
STATUS = "SUCCESS_POSTHOC_EUROPE_BIASED_INVENTORY_SENSITIVITY_NOT_SPATIAL_NULL"
EXPECTED_TOTAL = {
    "focal_butterflies":239,
    "base_added":14553,
    "augmented_added":14967,
    "reclassified_native":1403,
    "previously_unadded_new":1817,
    "affected_butterflies":81,
}
MARKERS = ("Lycaena dispar", "Minois dryas", "Thymelicus lineola")


def inspect_source(data: dict) -> None:
    if data.get("schema") != SCHEMA or data.get("status") != STATUS:
        raise ValueError("unrecognized source receipt")
    if data["source"].get("candidate_additional_host_links") != 1027:
        raise ValueError("candidate-host count mismatch")
    if data["source"].get("butterflies_with_candidate_hosts") != 81:
        raise ValueError("affected-species count mismatch")
    if data["source"].get("bce_accepted_edges_sha256") != (
        "f46d18703fe5cc8d97699af6f084f652a704e46b0f9e8c68e32dd70decc32c53"
    ):
        raise ValueError("BCE source crosswalk was not verified")
    b, a, diff = (data[k] for k in (
        "base", "BCE_augmented_partial", "change"
    ))
    expected = (
        (b["focal_butterflies"], a["focal_butterflies"]),
        (b["distinct_added_cells"], a["distinct_added_cells"]),
        (diff["old_added_cells_reclassified_native"],
         diff["previously_unadded_regions_newly_added"]),
    )
    if expected != ((239,239),(14553,14967),(1403,1817)):
        raise ValueError("original or augmented geographic envelope count drift")
    if diff["delta_U"] != 1817-1403:
        raise ValueError("source scenario net-added accounting failed")
    if len(data.get("affected_butterflies_sorted",[]))!=81:
        raise ValueError("paired species missing")


def species_turnover(row: dict) -> dict:
    old = int(row["base_added"])
    new = int(row["augmented_added"])
    removed = int(row["old_added_regions_reclassified_as_native"])
    added = int(row["previously_unadded_regions_newly_added"])
    retained = old - removed
    union = old + added
    if not (0 <= removed <= old and old-removed+added == new):
        raise ValueError(f"cannot reconcile per-butterfly added region identity: {row['species']}")
    out = {
        "species": row["species"],
        "BCE_candidate_host_links_added": int(row["BCE_candidate_host_links_added"]),
        "base_added": old,
        "augmented_added": new,
        "net_change": new-old,
        "retained_new_resource_cells": retained,
        "reclassified_native": removed,
        "previously_unadded_new": added,
        "union_new_resource_cells": union,
        "jaccard": retained / union if union else None,
        "is_turnover_despite_unchanged_count": bool(old == new and removed > 0 and added > 0),
        "both_reclassification_and_discovery": bool(removed > 0 and added > 0),
    }
    return out


def summarize(source:dict) -> dict:
    inspect_source(source)
    species = [species_turnover(r) for r in source["affected_butterflies_sorted"]]
    if len({x["species"] for x in species}) != 81:
        raise ValueError("duplicate augmented butterfly")
    add = sum(x["previously_unadded_new"] for x in species)
    rem = sum(x["reclassified_native"] for x in species)
    retained = sum(x["retained_new_resource_cells"] for x in species)
    union = sum(x["union_new_resource_cells"] for x in species)
    old_n = sum(x["base_added"] for x in species)
    new_n = sum(x["augmented_added"] for x in species)
    if (old_n, new_n, retained, union, rem, add) != (
        4628,5042,3225,6445,1403,1817
    ):
        raise ValueError("paired region identity drift among 81 augmented butterflies")
    valid = [x["jaccard"] for x in species if x["jaccard"] is not None]
    if len(valid)!=77:
        raise ValueError("nonempty species-level region union mismatch")
    positive = [x for x in species if x["net_change"]>0]
    negative = [x for x in species if x["net_change"]<0]
    zero = [x for x in species if x["net_change"]==0]
    double = [x for x in species if x["both_reclassification_and_discovery"]]
    if (len(positive),len(negative),len(zero),len(double)) != (36,35,10,54):
        raise ValueError("turnover direction distribution mismatch")
    hidden = [x for x in species if x["is_turnover_despite_unchanged_count"]]
    if [x["species"] for x in hidden] != ["Thymelicus lineola"]:
        raise ValueError("equal-count-but-shifted-identity case mismatch")
    extremes = {
        "most_positive_net":max(species, key=lambda x:(x["net_change"],x["species"])),
        "most_negative_net":min(species, key=lambda x:(x["net_change"],x["species"])),
        "equal_count_hidden_turnover":hidden[0],
    }
    if (
        extremes["most_positive_net"]["species"]!="Lycaena dispar"
        or extremes["most_negative_net"]["species"]!="Minois dryas"
    ):
        raise ValueError("extremal butterfly examples changed")
    full_retained = source["base"]["distinct_added_cells"]-rem
    full_union = source["base"]["distinct_added_cells"]+add
    if (full_retained,full_union)!=(13150,16370):
        raise ValueError("full-panel source contrast drift")
    return {
        "schema":"chocho_butterfly_bce_added_region_identity_turnover_v0.1",
        "status":"SUCCESS_DESCRIPTIVE_EUROPE_BIASED_HOST_INVENTORY_TURNOVER",
        "source":{
            "input_schema":source["schema"],
            "BCE_source_workflow_run":37795847333,
            "BCE_sensitivity_workflow_run":38988757959,
            "BCE_candidate_host_links_added":1027,
            "BCE_butterflies_modified":81,
            "focal_butterflies":239,
            "identity":"butterfly species x WGSRPD3 resource region, with first host support counted only once",
        },
        "affected_81":{
            "original_new_resource_cells":old_n,
            "augmented_new_resource_cells":new_n,
            "retained_in_both":retained,
            "reclassified_as_native":rem,
            "newly_added_by_augmentation":add,
            "union_of_either_source":union,
            "intersection_over_union_jaccard":retained/union,
            "discordant_fraction_of_union":(rem+add)/union,
            "nonempty_union_species":len(valid),
            "median_species_jaccard_nonempty_union":median(valid),
            "species_increase_net":len(positive),
            "species_decrease_net":len(negative),
            "species_equal_net":len(zero),
            "species_with_both_discovery_and_reclassification":len(double),
            "source_identity_discordance_overall":"This is item-specific overlap, NOT a p value or evidence of regional larval use.",
        },
        "full_239":{
            "original_new_resource_cells":14553,
            "augmented_new_resource_cells":14967,
            "retained_in_both":full_retained,
            "union_of_either_source":full_union,
            "intersection_over_union_jaccard":full_retained/full_union,
            "unchanged_host_inventory_species":158,
            "warning":"All 158 species without added BCE hosts are identical in the two scenarios, raising apparent full-panel overlap.",
        },
        "auditable_extremes":extremes,
        "species":species,
        "science_boundary":[
            "Source variation is Europe-biased and unidirectional. The 81 augmented butterflies are not an independently sampled global cohort.",
            "Sums compare geographic classifications under two catalogues and cannot show temporal ecological turnover.",
            "The Jaccard value is exact for the two mapped source scenarios but not an estimate of biological truth.",
            "Native and introduced WCVP statuses describe contemporary range-origin class, not dates or ability of larvae to complete development locally.",
            "Large source-induced shifts can occur with no total-count change, as demonstrated by Thymelicus lineola; no claim of real butterfly redistribution is permitted.",
            "Original HOSTS local spatial null was not recalibrated under extra BCE accepted host links.",
        ],
    }


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-json",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--output-species-csv",type=Path,required=True)
    args=ap.parse_args()
    src=json.loads(args.input_json.read_text(encoding="utf-8"))
    out=summarize(src)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    args.output_species_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_species_csv.open("w",encoding="utf-8",newline="") as h:
        rows=out["species"]
        w=csv.DictWriter(h,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    compact={k:v for k,v in out.items() if k!="species"}
    print(json.dumps(compact,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

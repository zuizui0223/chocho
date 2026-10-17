#!/usr/bin/env python3
"""Self-contained, standard-library butterfly larval-resource reanalysis.

This program is intentionally author-anonymous and uses only local inputs.
No network requests or personal identification are required. Reruns core
239-butterfly resource accounting, evidence-filtered 81-species host-list
supplementation, and exact butterfly x WGSRPD3 source-classification overlap.
MCMC estimands are supplied as source-qualified receipts, not silently rerun.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"code"))

from analyze_butterfly_shared_host_structural_null import build_model
from analyze_butterfly_climate_matched_plant_range_null import evaluate
from analyze_bce_resource_novelty_source_sensitivity import (
    load_candidate_host_links, host_geography, rows, analyze
)
from summarize_bce_added_resource_identity_turnover import summarize


def sha(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_inputs() -> dict:
    manifest=json.loads((ROOT/"SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("schema")!="butterfly_anonymous_offline_resource_package_v0.1":
        raise ValueError("not the accepted anonymised source manifest")
    files=manifest.get("sha256",{})
    if not files:
        raise ValueError("no verified files in archive")
    for relative,expected in files.items():
        if relative.startswith("/") or ".." in Path(relative).parts:
            raise ValueError("unsafe manifest entry")
        source=ROOT/relative
        if not source.is_file() or sha(source)!=expected:
            raise ValueError(f"source-file checksum differs: {relative}")
    return manifest


def run() -> dict:
    manifest=checked_inputs()
    data=ROOT/"inputs"
    d=json.loads((data/"level1_lookup.json").read_text(encoding="utf-8"))["units"]
    def csv_rows(relative):
        return rows(data/relative)

    host_rows=csv_rows("original/insect_host_accepted.csv")
    native_rows=csv_rows("original/native_extant_nondoubtful_wgsrpd3.csv")
    cont_rows=csv_rows("original/extant_nondoubtful_wgsrpd3.csv")
    desc=csv_rows("s1_resource_descriptors.csv")
    model=build_model(desc,host_rows,native_rows,cont_rows,d)
    observed=evaluate(model,model["introduced_sets"])
    for key,val in {
        "introduced_host_butterfly_region_incidences":75346,
        "introduced_incidences_in_native_supported_cells":46387,
        "redundant_introduced_host_incidences_in_new_cells":14406,
        "added_units":14553,
    }.items():
        if observed[key] != val:
            raise ValueError(f"Original HOSTS resource accounting changed: {key}")

    added,prov=load_candidate_host_links(data/"bce/accepted_host_edges.json")
    source=analyze(
        added,prov,
        csv_rows("anthropogenic_species_metrics.csv"),
        host_rows, host_geography(data/"original/native_extant_nondoubtful_wgsrpd3.csv"),
        host_geography(data/"original/extant_nondoubtful_wgsrpd3.csv"),
        host_geography(data/"bce/additional_native_wgsrpd3.csv"),
        host_geography(data/"bce/additional_contemporary_wgsrpd3.csv"),
    )
    turnover=summarize(source)
    if (
        source["base"]["native_absorbed_incidences"]!=46387
        or source["BCE_augmented_partial"]["candidate_incidences"]!=103349
        or source["BCE_augmented_partial"]["native_absorbed_incidences"]!=67740
        or source["BCE_augmented_partial"]["distinct_added_cells"]!=14967
        or turnover["affected_81"]["retained_in_both"]!=3225
        or turnover["affected_81"]["union_of_either_source"]!=6445
    ):
        raise ValueError("source sensitivity does not match audited study")
    # The published source crosswalk is validly evidence-filtered but Europe-
    # focused, so never elevate these source effects to worldwide biological
    # frequencies, and never inherit an original-HOSTS null p-value.
    original=source["base"]
    augmented=source["BCE_augmented_partial"]
    affected=turnover["affected_81"]
    sanitized={
        "schema":"anonymous_butterfly_resource_novelty_reproduced_v0.1",
        "status":"PASS_OFFLINE_EXACT_SOURCE_RECONSTRUCTION",
        "source_reconstruction_scope":"239 recorded butterflies; 365 WGSRPD3 resource regions",
        "original_HOSTS":{
            "candidate_incidences_H":original["candidate_incidences"],
            "native_absorbed_A":original["native_absorbed_incidences"],
            "novel_region_extra_support_R":original["redundant_new_cell_incidences"],
            "distinct_new_regions_U":original["distinct_added_cells"],
            "native_regions":original["native_cells"],
            "contemporary_regions":original["contemporary_cells"],
            "native_absorption_fraction":original["native_absorption_fraction"],
        },
        "BCE_partial_augmented":{
            "candidate_incidences_H":augmented["candidate_incidences"],
            "native_absorbed_A":augmented["native_absorbed_incidences"],
            "novel_region_extra_support_R":augmented["redundant_new_cell_incidences"],
            "distinct_new_regions_U":augmented["distinct_added_cells"],
            "native_regions":augmented["native_cells"],
            "contemporary_regions":augmented["contemporary_cells"],
            "native_absorption_fraction":augmented["native_absorption_fraction"],
        },
        "inventory_difference":{
            "native_overlap_share_delta_percentage_points":
                source["change"]["delta_native_absorption_fraction_percentage_points"],
            "old_new_regions_reclassified_native":
                source["change"]["old_added_cells_reclassified_native"],
            "unobserved_new_regions_identified":
                source["change"]["previously_unadded_regions_newly_added"],
            "net_added_region_difference":source["change"]["delta_U"],
        },
        "81_augmented_butterflies":{
            "source_classified_original_new_region_cells":affected["original_new_resource_cells"],
            "source_classified_augmented_new_region_cells":affected["augmented_new_resource_cells"],
            "retained_in_both":affected["retained_in_both"],
            "union_of_either":affected["union_of_either_source"],
            "jaccard":affected["intersection_over_union_jaccard"],
            "species_count_up":affected["species_increase_net"],
            "species_count_down":affected["species_decrease_net"],
            "species_count_same":affected["species_equal_net"],
        },
        "qualification":"Only original-HOSTS climate/geographic MCMC receipts supplied; not rerun after host augmentation.",
        "warnings":[
            "Catalogues differ, not historical butterfly ranges.",
            "This partially augmented host inventory is European, not a second globally representative sample.",
            "Potential larval resources are not local host feeding or larval fitness.",
            "The native-overlap partition and geographic novelty are not causal estimates.",
        ],
        "source_member_count":len(manifest["sha256"]),
    }
    outputs=ROOT/"reproduced"
    outputs.mkdir(exist_ok=True)
    (outputs/"reconstructed_resource_novelty.json").write_text(
        json.dumps(sanitized,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    with (outputs/"source_turnover_81_butterflies.csv").open("w",newline="",encoding="utf-8") as f:
        names=list(turnover["species"][0])
        w=csv.DictWriter(f,fieldnames=names)
        w.writeheader()
        w.writerows(turnover["species"])
    print(json.dumps(sanitized,ensure_ascii=False,indent=2))
    return sanitized


if __name__=="__main__":
    run()

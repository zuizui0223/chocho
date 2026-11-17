#!/usr/bin/env python3
"""Strict double-unseen botanical islands: true uniform sea CIRCUIT vs baseline.

Continues the same 104 WCVP island regions and 5,355 original accepted native
and introduced pair responses. Identical outcome-blind five island folds of
the earlier SOURCE-LOCKED shortest-water strict transfer are reused here;
ONLY the sea scalar is replaced with an independently computed undirected
effective resistance of the ORIGINAL GSHHG uniform finite-volume network.
This is a post-hoc transfer diagnostic, not a confirmatory introduction test.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as common
import audit_WCVP126_native_alien_double_unseen_islands_20261117 as strict

EXPECTED_CIRCUIT_PAIRS=5356
EXPECTED_MODEL_PAIRS=5355
EXPECTED_ISLANDS=104
EXPECTED_FOLDS=5

def csvrows(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("plants","shortest-water","uniform-circuit","uniform-circuit-receipt",
                "island-area","out-dir"):
        p.add_argument("--"+key,required=True,type=Path)
    args=p.parse_args()
    receipt=json.loads(args.uniform_circuit_receipt.read_text())
    if receipt["status"]!="SOURCE_VERIFIED_ORIGINAL_UNIFORM_MARINE_CIRCUIT_104_REGIONS_NOT_OCEAN_CURRENTS":
        raise ValueError("Original source graph provenance status changed")
    if receipt["source_original_104_unordered_pair_resistances"]!=EXPECTED_CIRCUIT_PAIRS:
        raise ValueError("Original source graph pair support changed")
    inputs=type("FrozenSources",(),{"plants":args.plants,"geometry":args.shortest_water})()
    pairs,not_defined=common.input_pair(inputs)
    if len(pairs)!=EXPECTED_MODEL_PAIRS or len(not_defined)!=1:
        raise ValueError("Original botanical pair response population changed")
    source=csvrows(args.uniform_circuit)
    if len(source)!=EXPECTED_CIRCUIT_PAIRS:
        raise ValueError("Original 5356 independently coast-anchored resistances changed")
    mapping={}
    for r in source:
        key=tuple(sorted((r["source_region"],r["destination_region"])))
        if key in mapping:
            raise ValueError("Duplicate undirected marine circuit region pair")
        v=float(r["log1p_normalized_uniform_circuit"])
        if not math.isfinite(v) or v<=0:
            raise ValueError("Nonfinite or nonpositive source circuit resistance")
        mapping[key]=v
    for x in pairs:
        pair=tuple(sorted((x["island1"],x["island2"])))
        if pair not in mapping:
            raise ValueError("Missing uniform resistance for WCVP botanical response pair")
        x["sea_detour"]=mapping[pair]
    area_rows=csvrows(args.island_area)
    area={x["code"]:float(x["area_km2"]) for x in area_rows}
    if len(area_rows)!=126 or len(area)!=126 or any(x<=0 for x in area.values()):
        raise ValueError("Original 126 regional source areas changed")
    native=strict.evaluate(pairs,"native",area)
    introduced=strict.evaluate(pairs,"alien",area)
    if native["n_double_unseen_test_pairs"]!=introduced["n_double_unseen_test_pairs"]:
        raise AssertionError("Not exactly same paired native/introduced double-unseen sample")
    if native["n_source_islands"]!=EXPECTED_ISLANDS or introduced["n_source_islands"]!=EXPECTED_ISLANDS:
        raise AssertionError("Original 104 island geography changed")
    if len(native["folds"])!=EXPECTED_FOLDS or len(introduced["folds"])!=EXPECTED_FOLDS:
        raise AssertionError("Strict fold count changed")
    for a,b in zip(native["folds"],introduced["folds"]):
        if a["held_island_codes"]!=b["held_island_codes"] or a["n_test_double_unseen_pairs"]!=b["n_test_double_unseen_pairs"]:
            raise AssertionError("Different native and introduced holdout islands")
    def extract(k):
        return native["macro_equal_fold_relative_mse_improvement_percent"][k], introduced["macro_equal_fold_relative_mse_improvement_percent"][k]
    report={
        "schema":"chocho_true_uniform_circuit_104islands_WCVP_native_alien_5fold_double_unseen_v1",
        "status":"SOURCE_VERIFIED_STRICT_BOTH_ISLANDS_UNSEEN_UNIFORM_OCEAN_CIRCUIT_NOT_CURRENT_OR_COLONISATION",
        "original_floral_source_run":44010011288,
        "original_circuit_source_run":44013658955,
        "original_one_island_holdout_run":44013825017,
        "original_shortest_water_double_unseen_run":44011411479,
        "graph_source_sha256":receipt["original_GSHHG_source_graph_SHA256"],
        "common_104_islands":EXPECTED_ISLANDS,
        "common_5355_flora_pairs":EXPECTED_MODEL_PAIRS,
        "strict_test_unique_double_unseen_pairs":native["n_double_unseen_test_pairs"],
        "fold_predefinition":"Identical outcome-blind SHA256 chocho-WCVP126-double-unseen-v1|region hash, five bins, no response-selected fold assignment",
        "uniform_circuit_added_vs_base_pct":{"native":extract("G_vs_B_pct")[0],"introduced":extract("G_vs_B_pct")[1]},
        "contemporary_port_count_added_vs_base_pct":{"native":extract("P_vs_B_pct")[0],"introduced":extract("P_vs_B_pct")[1]},
        "native":native,"introduced":introduced,
        "limits":[
            "This is a post-hoc same-WCVP-data transfer stress, not independent discovery replication.",
            "Both endpoints of each tested pair are absent from training but their contemporary WCVP richness marginal counts enter the predictors.",
            "Only 5 held-island-fold groups; descriptive bootstrap over folds cannot yield independent-group confirmatory intervals.",
            "The circuit is GSHHG undirected finite-volume electrical conductance, not true directional ocean current or observed transport.",
            "Original WCVP taxa are current native or introduced listings, not observed historical arrivals.",
            "Contemporary WPI port count product is not verified vessel OD.",
        ],
    }
    args.out_dir.mkdir(parents=True,exist_ok=True)
    out=args.out_dir/"WCVP126_native_alien_uniform_circuit_BOTH_islands_unseen_5fold_20261117.json"
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print("WCVP126_CIRCUIT_BOTH_UNSEEN_SOURCE_LOCKED",json.dumps({
         "source_islands":EXPECTED_ISLANDS,"common_pairs":EXPECTED_MODEL_PAIRS,
         "strict_test_pairs":report["strict_test_unique_double_unseen_pairs"],
         "uniform_circuit_effects_pct":report["uniform_circuit_added_vs_base_pct"],
         "native_descriptive5fold95":native["descriptive_5fold_bootstrap95_not_independent_validation"]["G_vs_B_pct"],
         "alien_descriptive5fold95":introduced["descriptive_5fold_bootstrap95_not_independent_validation"]["G_vs_B_pct"],
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()

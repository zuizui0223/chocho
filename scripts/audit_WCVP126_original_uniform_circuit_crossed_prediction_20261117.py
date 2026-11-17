#!/usr/bin/env python3
"""Exact same 5355 accepted WCVP native/alien island pairs: uniform-CIRCUIT vs ports.

The source GSHHG circuit is solved from original source conductance and
104 verified nonport coastal terminal distributions. This replays the
original native/alien botanical leave-one-island-out scoring framework with
identical same original botanical pair, response, richness/area, port index,
training-only standardization and ridge alpha20, replacing ONLY the former
sea shortest-water detour with the source-equivalent full uniform-ocean
electrical effective resistance.

No claim of currents, ship/seed movements, changes in butterfly range.
"""
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path

import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as original

ORIGINAL_CIRCUIT_CI_SOURCE="pending exact uniform source"
ORIGINAL_GSHHG_SHORTEST_MARINE_SOURCE_CI=44010673580
ORIGINAL_NOAA_2005_SOURCE_CI=44012093616
ORIGINAL_COMMON_ISLANDS=104
ORIGINAL_COMMON_PAIRS=5355
ORIGINAL_ELECTRIC_PAIRS=5356

def read(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser()
    for k in ("plants","shortest-water","uniform-circuit","uniform-circuit-receipt",
              "original-sea-audit","original-shipping-audit","island-area","out-dir"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    circ=json.loads(a.uniform_circuit_receipt.read_text())
    if (circ["status"]!="SOURCE_VERIFIED_ORIGINAL_UNIFORM_MARINE_CIRCUIT_104_REGIONS_NOT_OCEAN_CURRENTS" or
        circ["source_original_104_unordered_pair_resistances"]!=ORIGINAL_ELECTRIC_PAIRS or
        circ["source_component_marine_nodes"]!=189550):
        raise RuntimeError("Original GSHHG finite-volume uniform circuit source identity changed")
    original_inputs=type("FrozenSourcePaths",(object,),{
        "plants":a.plants,"geometry":a.shortest_water,
    })()
    pairs,unknown=original.input_pair(original_inputs)
    if len(pairs)!=ORIGINAL_COMMON_PAIRS or len(unknown)!=1:
        raise RuntimeError(f"Original same-sample 104×104 accepted WCVP native/alien Jaccard source changed: {len(pairs)}")
    source_circuit=read(a.uniform_circuit)
    if len(source_circuit)!=ORIGINAL_ELECTRIC_PAIRS:
        raise RuntimeError("Source original circuit 5356 pair result changed")
    R={}
    for r in source_circuit:
        key=tuple(sorted((r["source_region"],r["destination_region"])))
        if key in R:raise RuntimeError("Duplicate source circuit region pair")
        val=float(r["log1p_normalized_uniform_circuit"])
        if not math.isfinite(val) or val<=0:
            raise RuntimeError("Source uniform circuit normalized scalar invalid")
        R[key]=val
    oldgeo=[e["sea_detour"] for e in pairs]
    for e in pairs:
        key=tuple(sorted((e["island1"],e["island2"])))
        if key not in R:raise RuntimeError("Source source-circuit/source original pair mismatch")
        e["sea_detour"]=R[key]
    b=read(a.island_area)
    areas={r["code"]:float(r["area_km2"]) for r in b}
    if len(areas)!=126 or len(b)!=126:
        raise RuntimeError("Original WCVP126 island area census changed")
    native=original.evaluate(pairs,"native",areas)
    alien=original.evaluate(pairs,"alien",areas)
    desc={
      "source_count_104_islands":ORIGINAL_COMMON_ISLANDS,
      "source_common_native_alien_valid_pairs":ORIGINAL_COMMON_PAIRS,
      "source_circuit_eligible_region_pairs":ORIGINAL_ELECTRIC_PAIRS,
      "native_B_to_G_circuit_improvement_pct":
           native["relative_delta_MSE_percent_positive_is_improvement"]["G_vs_B_pct"],
      "alien_B_to_G_circuit_improvement_pct":
           alien["relative_delta_MSE_percent_positive_is_improvement"]["G_vs_B_pct"],
      "native_after_port_add_circuit_pct":
           native["relative_delta_MSE_percent_positive_is_improvement"]["G_added_after_P_pct"],
      "alien_after_port_add_circuit_pct":
           alien["relative_delta_MSE_percent_positive_is_improvement"]["G_added_after_P_pct"],
      "native_port_added_after_circuit_pct":
           native["relative_delta_MSE_percent_positive_is_improvement"]["P_added_after_G_pct"],
      "alien_port_added_after_circuit_pct":
           alien["relative_delta_MSE_percent_positive_is_improvement"]["P_added_after_G_pct"]
    }
    doc={
      "schema":"chocho_104island_5355_native_alien_true_uniform_circuit_VS_modern_port_proxy_same_holdout_v1",
      "status":"SOURCE_LOCKED_BIOTIC_PREDICTION_WITH_TRUE_UNIFORM_GSHHG_CIRCUIT_CONTROL_NOT_OCEAN_CURRENTS",
      "original_WCVP_source_CI":44010011288,
      "original_GSHHG_coastline_source_CI":ORIGINAL_GSHHG_SHORTEST_MARINE_SOURCE_CI,
      "original_2020_AIS_source_CI":44010307551,
      "original_uniform_finite_volume_graph_sha256":
           circ["original_GSHHG_source_graph_SHA256"],
      "source_original_uniform_finvol_104_island_circuit_solver_receipt":circ,
      "model_variable_natural_geography":
          "log1p(actual source-uniform GSHHG finite-volume electrical effective resistance / median of all original 104 same-source resistances)",
      "model_variable_anthropogenic_port":
          "product log1p(current 2023 WPI port counts) of both exact botanical island regions; NOT shipping voyages",
      "original_baseline_same_WCVP_pairs":
          "Geodesic log distance plus quadratic, TDWG Level1 equality, latitude difference, native/alien accepted source species richness min/max and island area min/max",
      "original_held_source_island_dyad_folds":ORIGINAL_COMMON_ISLANDS,
      "source_original_equal_held_island_validation_ridge_alpha":20,
      "native_flora":native,
      "introduced_flora":alien,
      "contrasts":desc,
      "limitations":[
        "Full source-defined uniform-ocean circuit resistance is NOT GSHHG shortest-water path and NOT a physically directional model of NOAA currents.",
        "Nonport shoreline sampled 1–3 times per TDWG region, not a true seed spawning site; archipelago topology can be coarsened.",
        "Port count proxy is NOT observed vessel OD, especially 2020 strict WPI AIS links only 7/7875 full source botanical island pairs.",
        "The current state native/alien accepted WCVP floras are co-reported and not historically observed seed arrival.",
        "This test controls distance and richness for co-recorded modern plant distributions; ecological colonization/cargo history and butterflies remain unobserved.",
        "The same 5355 existing WCVP native/alien pair outcomes are used in both shortest-water and electrical tests; results are not independent biological data replications.",
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"WCVP126_native_alien_true_uniform_circuit_vs_port_held_island_validation.json").write_text(
       json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("WCVP126_TRUE_UNIFORM_CIRCUIT_CROSSED_PREDICTION",json.dumps({
       "status":doc["status"],
       "summary":desc,
       "native95_G":native["descriptive_island_bootstrap_95pct_MSE_delta"]["G_vs_B_pct"],
       "alien95_G":alien["descriptive_island_bootstrap_95pct_MSE_delta"]["G_vs_B_pct"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()

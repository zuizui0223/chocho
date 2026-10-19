#!/usr/bin/env python3
"""Finer plant-climate/geographic conditional null for butterfly host resources.

Data: source-frozen original HOSTS/WCVP/CHELSA 239 butterflies + 1027
partial BCE/Clarke host links (81 butterfly taxa, 679 additional accepted
plant IDs). All source geographic inputs remain author-only.

The only observed response is potential newly mapped larval host-resource
climate analogue fraction (O/A); NOT insect occupancy, colonisation, health,
larval survival, population growth, selection, or conservation benefit.

See provenance and stopping gates:
docs/exploratory/BCE_FINER_CLIMATE_SPATIAL_FILTER_PROTOCOL_20261019.md

Only pure summary/complete numeric MCMC draws are saved; no WCVP plant
distribution matrices, butterfly accepted-host edges or TDWG geometry are
published in the result.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_bce_augmented_plant_network_climate_null_20261019 import (
    ALL_RESPONSE_KEYS, EXPECTED_INTRODUCED_PLANT_REGION_EDGES,
    SEEDS, marginal_signature, source_augmented_model, analog_mask, metrics,
)
from analyze_climate_analog_fixedmargin_null_20261019 import ess_approx

MODES=(
    "plant_bio1_bio12_3x3",
    "plant_bio1_bio7_bio12_bio15_2x2x2x2",
    "plant_bio1_bio12_2x2_local1500km",
)
MODEL_MODES={
    "plant_bio1_bio12_3x3":"temperature_precipitation_3x3",
    "plant_bio1_bio7_bio12_bio15_2x2x2x2":"four_bioclim_variables_2x2x2x2",
    "plant_bio1_bio12_2x2_local1500km":"temperature_precipitation_2x2",
}
N_DRAWS=199
BURN_FACTOR=10
SPACING_FACTOR=4
CALIPER_KM=1500
EXPECTED_OBS={
    "0.05":(9949,14961),
    "0.1":(7370,14961),
    "matched_0.05":(9903,14903),
    "matched_0.1":(7331,14903),
}
MIN_VIABLE_EDGE_FRACTION=0.75
MIN_ORIGINAL_EDGES_MOVED_CLIMATE=0.10
MIN_ORIGINAL_EDGES_MOVED_LOCAL=0.05
MIN_ACCEPTED_PER_SPACING=100
MIN_ESS=30


def source_stratum_setup(source_root,model,regions,climate,original,mode):
    """Original versioned climate and geometry functions; no hand-coded bins."""
    sys.path.insert(0,str(source_root/"code"))
    from analyze_butterfly_climate_matched_plant_range_null import (
        climate_strata, grouped_introduced_edges,
    )
    classes,metadata=climate_strata(
        climate,regions,model["level1"],MODEL_MODES[mode])
    groups,cumulative,structure=grouped_introduced_edges(
        original,model["host_families"],model["level1"],classes)
    total=structure["total_introduced_plant_region_edges"]
    assert total==EXPECTED_INTRODUCED_PLANT_REGION_EDGES
    assert structure["climate_covered_edges"]==total
    assert structure["viable_fraction_of_total"]>=MIN_VIABLE_EDGE_FRACTION
    spatial={}
    if mode=="plant_bio1_bio12_2x2_local1500km":
        from analyze_butterfly_spatially_local_plant_range_null import (
            load_geography, neighbour_index, origin_tether,
        )
        coords,geometry_meta=load_geography(
            source_root/"inputs/level3.geojson",regions,model["level1"])
        near=neighbour_index(coords,CALIPER_KM)
        tether=origin_tether(original,near)
        spatial={
            "neighbours":near,"tether":tether,
            "spatial_provenance":{
                **geometry_meta,
                "radius_km":CALIPER_KM,
                "mean_geographic_neighbor_count":
                    sum(len(v) for v in near.values())/len(near),
                "min_neighbor_count":min(len(v) for v in near.values()),
                "original_introduced_plant_footprint_tethered":True,
                "distance_is_not_insect_dispersal_capacity":True,
            }
        }
    return classes,groups,cumulative,structure,metadata,spatial


def run_chain(source_root,model,regions,climate,observed_masks,original,mode,seed,
              *,draws=N_DRAWS,burn=BURN_FACTOR,spacing=SPACING_FACTOR):
    from analyze_butterfly_shared_host_structural_null import attempt_range_swaps
    from analyze_butterfly_spatially_local_plant_range_null import local_swap
    plants={host:set(cells) for host,cells in original.items()}
    classes,groups,cumulative,structure,metadata,spatial=source_stratum_setup(
        source_root,model,regions,climate,plants,mode)
    E=sum(len(cells) for cells in plants.values())
    assert E==EXPECTED_INTRODUCED_PLANT_REGION_EDGES
    baseline_signature=marginal_signature(model,original,classes)
    assert marginal_signature(model,plants,classes)==baseline_signature
    rng=random.Random(seed)

    def do_swaps(n):
        if spatial:
            return local_swap(
                groups,cumulative,plants,model["native_sets"],
                spatial["neighbours"],spatial["tether"],n,rng)
        return attempt_range_swaps(
            groups,cumulative,model["native_sets"],plants,n,rng)

    accepted_burn=do_swaps(burn*E)
    moved=sum(len(plants[h]^original[h]) for h in original)/(2*E)
    assert marginal_signature(model,plants,classes)==baseline_signature
    history={key:[] for key in ALL_RESPONSE_KEYS}
    accepted_per_draw=[]
    for draw in range(draws):
        accepted=do_swaps(spacing*E)
        accepted_per_draw.append(accepted)
        all_out,all_added=metrics(model,plants,observed_masks)
        assert set(all_out)==set(ALL_RESPONSE_KEYS)
        for key,row in all_out.items():
            history[key].append({
                "aligned":row["aligned"],
                "evaluable_added":row["evaluable_added"],
                "aligned_fraction":row["aligned_fraction"],
                "all_focal_butterflies_added_regions":all_added,
            })
        if draw in (0,draws//2,draws-1):
            if marginal_signature(model,plants,classes)!=baseline_signature:
                raise ValueError("PLANT_FAMILY_OR_CLIMATE_OR_SPATIAL_MARGIN_DRIFT")
    summary={}
    for key,draws_records in history.items():
        values=[r["aligned_fraction"] for r in draws_records]
        denom=[r["evaluable_added"] for r in draws_records]
        summary[key]={
            "ESS_diagnostic":ess_approx(values),
            "null_median":float(np.median(values)),
            "null_fraction_95_percentile":[
                float(np.quantile(values,q)) for q in (.025,.975)],
            "null_denominator_median":float(np.median(denom)),
            "null_denominator_range":[min(denom),max(denom)],
        }
    return {
        "seed":seed,
        "draws":draws,
        "mode":mode,
        "plant_edges":E,
        "stratum_structure":structure,
        "climate_strata":metadata,
        "geography":spatial.get("spatial_provenance"),
        "burn_proposals":burn*E,
        "burn_accepted_swaps":accepted_burn,
        "fraction_original_plant_edges_moved_after_burn":moved,
        "spacing_proposals":spacing*E,
        "accepted_swaps_per_spacing_mean":float(np.mean(accepted_per_draw)),
        "summary":summary,
        "full_draw_history":history,
        "provenance":"Source exact and native exclusions and all host×stratum and region×plant-family degrees checked after burn and at first/mid/last draws"
    }


def quality_and_contrast(observed,chains,mode):
    assert len(chains)==2 and [c["seed"] for c in chains]==list(SEEDS)
    threshold=(MIN_ORIGINAL_EDGES_MOVED_LOCAL if "local1500km" in mode
               else MIN_ORIGINAL_EDGES_MOVED_CLIMATE)
    qg={};out={}
    for key in ALL_RESPONSE_KEYS:
        stats=[c["summary"][key] for c in chains]
        med=[x["null_median"] for x in stats]
        e=[x["ESS_diagnostic"]["ess_approx"] for x in stats]
        std=[x["ESS_diagnostic"]["std"] for x in stats]
        ess_ok=all(v>=MIN_ESS for v in e)
        moved_ok=all(c["fraction_original_plant_edges_moved_after_burn"]>=threshold
                     for c in chains)
        acceptance=all(c["accepted_swaps_per_spacing_mean"]>=MIN_ACCEPTED_PER_SPACING
                       for c in chains)
        concordance=abs(med[0]-med[1])<=2*max(std)
        structure_ok=all(
            c["stratum_structure"]["viable_fraction_of_total"]>=MIN_VIABLE_EDGE_FRACTION
            for c in chains)
        qualified=ess_ok and moved_ok and acceptance and concordance and structure_ok
        qg[key]={
            "qualified":qualified,
            "ESS_each_ge_30":ess_ok,
            "burn_moved_fraction_pass":moved_ok,
            "accepted_per_spacing_ge_100":acceptance,
            "chain_medians_agree":concordance,
            "viable_plant_edges_ge_75_percent":structure_ok,
            "min_burn_moved_fraction_gate":threshold,
        }
        obs=observed[key]["aligned_fraction"]
        full_tail=[
            (1+sum(v["aligned_fraction"]>=obs for v in c["full_draw_history"][key]))/
            (len(c["full_draw_history"][key])+1) for c in chains
        ] if qualified else None
        out[key]={
            "observed":observed[key],
            "median_null_fraction_by_chain":med,
            "95pct_null_fraction_by_chain":[x["null_fraction_95_percentile"] for x in stats],
            "residual_percentage_points_by_chain":[100*(obs-v) for v in med],
            "null_denominator_median_by_chain":[x["null_denominator_median"] for x in stats],
            "null_denominator_range_by_chain":[x["null_denominator_range"] for x in stats],
            "ESS_approx_by_chain":e,
            "right_tail_serial_MCMC_not_exact_randomization_p_by_chain":full_tail,
            "qualified":qualified,
        }
    return {
        "status":"QUALIFIED_EXPLORATORY_SOURCE_FIXED_FINE_CLIMATE_NULL"
                 if all(v["qualified"] for v in qg.values()) else
                 "NOT_EVALUABLE_MIXING_OR_GEOGRAPHIC_LOCALITY",
        "quality":qg,
        "contrasts":out,
        "chains":chains,
    }


def run(source_root,*,modes=MODES,draws=N_DRAWS,burn=BURN_FACTOR,spacing=SPACING_FACTOR):
    source_root=Path(source_root)
    model,regions,climate,provenance=source_augmented_model(source_root)
    assert len(model["species"])==239 and len(regions)==365
    masks,eligible=analog_mask(model,regions,climate)
    assert len(eligible)==235
    original={h:set(x) for h,x in model["introduced_sets"].items()}
    assert sum(map(len,original.values()))==EXPECTED_INTRODUCED_PLANT_REGION_EDGES
    observed,source_added=metrics(model,original,masks)
    assert source_added==14965
    for k,(O,A) in EXPECTED_OBS.items():
        assert (observed[k]["aligned"],observed[k]["evaluable_added"])==(O,A),k
    modes_out={}
    for mode in modes:
        assert mode in MODES
        chains=[
            run_chain(source_root,model,regions,climate,masks,original,mode,seed,
                      draws=draws,burn=burn,spacing=spacing)
            for seed in SEEDS
        ]
        modes_out[mode]=quality_and_contrast(observed,chains,mode)
    return {
        "schema":"chocho_BCE_239_butterflies_finer_plant_climate_spatial_MCMC_v0.1",
        "status":"QUALIFIED_EXPLORATORY_ALL_REQUESTED_MODES"
                 if all(v["status"]=="QUALIFIED_EXPLORATORY_SOURCE_FIXED_FINE_CLIMATE_NULL"
                        for v in modes_out.values()) else
                 "PARTIAL_OR_NOT_EVALUABLE",
        "source":"Exact author-only WCVP/HOSTS/BCE/CHELSA/TDWG permitted analysis",
        "source_integrity":provenance,
        "scope":{
            "butterflies":239,
            "climate_eligible_augmented_butterflies":235,
            "same_original_climate_eligible_butterflies":234,
            "additional_accepted_butterfly_host_links":1027,
            "added_accepted_plant_names":679,
            "source_introduced_plant_WGSRPD3_links":37869,
            "climate_eligible_new_resource_cells":14961,
            "unavailable_climate_new_regions_for_Coenonympha_oedippus":["ANT","HMD"],
        },
        "observed":observed,
        "seeds":list(SEEDS),
        "draws_per_chain":draws,
        "burn_attempts_per_plant_edge":burn,
        "spacing_attempts_per_plant_edge":spacing,
        "modes":modes_out,
        "claim_gates":[
            "Finer climate and geographic caliper nulls were selected after positive 2x2 residual was observed; exploratory not confirmatory preregistration.",
            "Area cells are unequal; plant climatic classes and representative centroids do not directly estimate plant or butterfly physiological niche.",
            "Spatial radius is not insect dispersal distance; original footprint tether defines a constrained reachable MCMC state space.",
            "Serial MCMC right tails are not independent permutation p-values.",
            "No observed butterfly occupancy, larval performance, climate adaptation, reproductive recruitment or causal effects of plant introductions.",
            "Host catalogue supplement is Europe-biased, no independent global host source replication.",
            "Full original accepted plant region source remains private for rights; only derived numeric MCMC chains are uploaded."
        ]
    }


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--source-root",type=Path,required=True)
    a.add_argument("--output-json",type=Path,required=True)
    a.add_argument("--mode",nargs="+",choices=MODES,default=list(MODES))
    a.add_argument("--draws",type=int,default=N_DRAWS)
    a.add_argument("--burn-factor",type=int,default=BURN_FACTOR)
    a.add_argument("--spacing-factor",type=int,default=SPACING_FACTOR)
    args=a.parse_args()
    result=run(args.source_root,modes=tuple(args.mode),draws=args.draws,
               burn=args.burn_factor,spacing=args.spacing_factor)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",
                                encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "source":"239 butterflies; BCE partial +1027 hosts; 37869 plant region links",
        "observed":result["observed"],
        "modes":{
            k:{
                "status":v["status"],
                "quality":v["quality"],
                "comparisons":v["contrasts"],
                "source_swap_viability":v["chains"][0]["stratum_structure"],
                "source_spatial_geography":v["chains"][0]["geography"],
                "moved_original_edge_frac":[c["fraction_original_plant_edges_moved_after_burn"]
                                             for c in v["chains"]],
                "mean_accepted_swaps":[c["accepted_swaps_per_spacing_mean"]
                                        for c in v["chains"]],
            } for k,v in result["modes"].items()
        }
    },sort_keys=True))

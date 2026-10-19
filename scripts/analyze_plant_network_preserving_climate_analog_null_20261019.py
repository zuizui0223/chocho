#!/usr/bin/env python3
"""239-butterfly climate analogue null preserving *plant-host network* structure.

The exact HOSTS/WCVP/CHELSA input archive is author-only. This program uses
existing source-model SHA gates, never uploads/redistributes plant-level inputs.

Two prospective nulls (protocol fixed before this full-chain run):
  1. introduced plant × region swaps within plant family × WGSRPD Level1;
  2. the same PLUS continent-local BIO1×BIO12 2×2 climate bins.

Each preserves butterfly-host edges, plant native ranges, host introduction
degree within allowed groups, and introduced host count by plant family×region.
Because distinct butterfly resource union sizes can change under swaps, the
primary response is the FRACTION of *evaluable* newly added butterfly×region
cells inside each butterfly's original native-host climatic rectangle, not
raw aligned-cell counts. Not a butterfly colonisation or demographic model.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import numpy as np

from analyze_butterfly_climate_analog_macro_20261019 import (
    bits_to_indices, qbounds, source_model,
)
from analyze_butterfly_shared_host_structural_null import (
    attempt_range_swaps, make_groups,
)
from analyze_butterfly_climate_matched_plant_range_null import (
    climate_strata, grouped_introduced_edges,
)
from analyze_climate_analog_fixedmargin_null_20261019 import ess_approx

SEEDS=(20261019,20261020)
QS=(0.05,0.10)
MODES=("plant_family_continent","plant_family_continent_climate_bio1_bio12_2x2")
BURN_FACTOR=10
SPACING_FACTOR=4
DRAW_COUNT=199
EXPECTED_OBS={"0.05":{"aligned":9451,"evaluable_added":14481},
              "0.1":{"aligned":6948,"evaluable_added":14481}}
EXPECTED_PLANT_EDGES=30373
MIN_VIABLE_PLANT_EDGE_SHARE=.75
MIN_ESS_PER_CHAIN=30
MIN_BURN_MOVED=.10
MIN_ACCEPTED_SPACING=50


def source_analog_masks(model,regions,climate,raw):
    # Only the original native HOSTS range builds each rectangle.
    # Never recompute the rectangle from the shuffled added resource.
    covered={i for i,code in enumerate(regions) if code in climate}
    assert len(covered)==365
    masks={q:{} for q in QS}
    supported=[]
    for i,row in enumerate(raw):
        native=set(bits_to_indices(row["native_bits"],len(regions)))
        native_climate=sorted(native & covered)
        if len(native_climate)<10:
            continue
        supported.append(i)
        for q in QS:
            lo1,hi1=qbounds([climate[regions[j]][0] for j in native_climate],q)
            lo12,hi12=qbounds([climate[regions[j]][2] for j in native_climate],q)
            b=0
            for j in covered:
                values=climate[regions[j]]
                if lo1<=values[0]<=hi1 and lo12<=values[2]<=hi12:
                    b |= 1<<j
            masks[q][i]=b
    assert len(supported)==234
    return masks, supported


def butterfly_metrics(model,plants,climate_masks):
    host_mask={host:sum(1<<j for j in indices)
               for host,indices in plants.items()}
    added_by_species=[]
    for i,hosts in enumerate(model["species_hosts"]):
        combined=0
        for host in hosts:
            combined |= host_mask.get(host,0)
        added_by_species.append(combined & ~model["native_species"][i])
    out={}
    for q,masks in climate_masks.items():
        gained=sum(added_by_species[i].bit_count() for i in masks)
        overlap=sum((added_by_species[i]&mask).bit_count()
                    for i,mask in masks.items())
        assert gained>0 and 0<=overlap<=gained
        out[str(q)]={
            "aligned":overlap,
            "evaluable_added":gained,
            "aligned_fraction":overlap/gained,
            "butterflies_with_any_added":sum(added_by_species[i]!=0 for i in masks),
        }
    return out, sum(v.bit_count() for v in added_by_species)


def signature(model,plants,climate_classes=None):
    # More demanding invariant assertions than the swap proposal alone:
    # plant × area(/climate) introduction degrees and
    # native-excluded introduced plant-family × region degrees.
    byplant=Counter()
    byregion=Counter()
    byfamily_area=Counter()
    count=0
    for host,cells in plants.items():
        fam=model["host_families"][host]
        forbidden=model["native_sets"].get(host,set())
        assert not set(cells)&forbidden, ("native-host-cell",host)
        for region in cells:
            area=model["level1"][region]
            geo=(area,climate_classes.get(region) if climate_classes else None)
            byplant[(host,geo)]+=1
            byregion[(region,fam)]+=1
            byfamily_area[(fam,geo)]+=1
            count+=1
    return byplant,byregion,byfamily_area,count


def compare_invariants(baseline,current):
    for i,(a,b) in enumerate(zip(baseline,current)):
        if a!=b:
            raise ValueError(f"PLANT_MARGIN_DRIFT dimension={i} example={list((a-b).items())[:4] if i<3 else (a,b)}")
    assert current[3]==EXPECTED_PLANT_EDGES


def make_mode_groups(model,climate,regions,plants,mode):
    E=sum(map(len,plants.values()))
    assert E==EXPECTED_PLANT_EDGES
    if mode=="plant_family_continent":
        groups,cumulative,total_edges=make_groups(
            model["native_sets"],plants,model["host_families"],model["level1"])
        assert total_edges==E
        classes={}
        metrics={
            "covered_fraction":1.0,
            "viable_edge_fraction":(cumulative[-1] if cumulative else 0)/E,
            "viable_groups":len(groups),
            "conditioning":"plant_family_x_continent"
        }
    elif mode=="plant_family_continent_climate_bio1_bio12_2x2":
        classes,meta=climate_strata(
            climate,regions,model["level1"],"temperature_precipitation_2x2")
        groups,cumulative,m=grouped_introduced_edges(
            plants,model["host_families"],model["level1"],classes)
        metrics={
            "covered_fraction":m["climate_covered_fraction"],
            "viable_edge_fraction":m["viable_fraction_of_total"],
            "viable_groups":m["viable_group_count"],
            "conditioning":"plant_family_x_continent_x_bio1_bio12_2x2",
            "non_swappable_groups":m["nonswappable_group_count"],
        }
    else:
        raise ValueError(mode)
    if metrics["covered_fraction"]<.95 or metrics["viable_edge_fraction"]<MIN_VIABLE_PLANT_EDGE_SHARE:
        raise ValueError("PLANT_CLIMATE_SUPPORT_OR_SWAPPABILITY_GATE")
    return groups,cumulative,classes,metrics


def run_chain(model,regions,climate,original,climate_masks,mode,seed):
    # Each seed starts from unmodified source; groups modify the same plant sets
    # as the proposal state. We NEVER seed chain 2 with chain 1.
    plants={k:set(v) for k,v in original.items()}
    groups,cumulative,classes,structure=make_mode_groups(
        model,climate,regions,plants,mode)
    baseline=signature(model,original,classes)
    compare_invariants(baseline,signature(model,plants,classes))
    E=baseline[3]
    rng=random.Random(seed)
    burn_acc=attempt_range_swaps(
        groups,cumulative,model["native_sets"],plants,BURN_FACTOR*E,rng)
    moved=sum(len(plants[h]^original[h]) for h in plants)/(2*E)
    compare_invariants(baseline,signature(model,plants,classes))
    samples={str(q):[] for q in QS}
    spacing_accepted=[]
    for i in range(DRAW_COUNT):
        accepted=attempt_range_swaps(
            groups,cumulative,model["native_sets"],plants,SPACING_FACTOR*E,rng)
        spacing_accepted.append(accepted)
        metrics,all_gained=butterfly_metrics(model,plants,climate_masks)
        for q in QS:
            samples[str(q)].append({
                "aligned":metrics[str(q)]["aligned"],
                "evaluable_added":metrics[str(q)]["evaluable_added"],
                "aligned_fraction":metrics[str(q)]["aligned_fraction"],
                "all_species_added":all_gained,
            })
        if i in (0,DRAW_COUNT-1,DRAW_COUNT//2):
            compare_invariants(baseline,signature(model,plants,classes))
    summaries={}
    for q in QS:
        k=str(q)
        draws=[x["aligned_fraction"] for x in samples[k]]
        synth=ess_approx(draws)
        srt=sorted(draws)
        summaries[k]={
            "n_draws":len(draws),
            "mean_fraction":float(np.mean(draws)),
            "median_fraction":float(median(draws)),
            "sd_fraction":float(np.std(draws,ddof=1)),
            "percentile_95_fraction":[float(np.quantile(draws,.025)),
                                      float(np.quantile(draws,.975))],
            "ess_diagnostic":synth,
            "mean_null_added_cells":float(np.mean(
                [x["evaluable_added"] for x in samples[k]])),
            "min_null_added_cells":min(x["evaluable_added"] for x in samples[k]),
            "max_null_added_cells":max(x["evaluable_added"] for x in samples[k]),
        }
    return {
        "seed":seed,
        "mode":mode,
        "structure":structure,
        "burn_attempts":BURN_FACTOR*E,
        "accepted_burn":burn_acc,
        "fraction_plant_edge_positions_changed_after_burn":moved,
        "accepted_swaps_per_spacing_mean":float(np.mean(spacing_accepted)),
        "spacing_attempts":SPACING_FACTOR*E,
        "summaries":summaries,
        "draws":samples,
        "integrity_checks":"Original plant×region family degrees and plant×stratum degrees checked after burn and 3 times during sampling, native exclusions verified.",
    }


def summarize_mode(chains,observed,mode):
    quality={}
    contrasts={}
    assert len(chains)==2
    for q in QS:
        k=str(q)
        a,b=(c["summaries"][k] for c in chains)
        moved=all(c["fraction_plant_edge_positions_changed_after_burn"]>=MIN_BURN_MOVED
                  for c in chains)
        ess=all(c["summaries"][k]["ess_diagnostic"]["ess_approx"]>=MIN_ESS_PER_CHAIN
                for c in chains)
        mixing=all(c["accepted_swaps_per_spacing_mean"]>=MIN_ACCEPTED_SPACING
                   for c in chains)
        median_agree=abs(a["median_fraction"]-b["median_fraction"])<=\
            2*max(a["sd_fraction"],b["sd_fraction"])
        viable=all(c["structure"]["viable_edge_fraction"]>=MIN_VIABLE_PLANT_EDGE_SHARE
                   for c in chains)
        ok=moved and ess and mixing and median_agree and viable
        quality[k]={
            "qualified":ok,"source_structure_viable":viable,
            "each_burn_moves_ge_0_10":moved,
            "each_ESS_ge_30":ess,
            "spacing_accepted_each_ge_50":mixing,
            "chain_null_median_agreement":median_agree,
        }
        obj=observed[k]
        obs_frac=obj["aligned_fraction"]
        right_tails=[
            (1+sum(d["aligned_fraction"]>=obs_frac for d in c["draws"][k]))/
            (len(c["draws"][k])+1)
            for c in chains
        ]
        contrasts[k]={
            "observed_evaluable_added":obj["evaluable_added"],
            "observed_aligned":obj["aligned"],
            "observed_aligned_fraction":obs_frac,
            "null_fraction_median_by_chain":[a["median_fraction"],b["median_fraction"]],
            "null_fraction_95pct_by_chain":[
                a["percentile_95_fraction"],b["percentile_95_fraction"]],
            "observed_minus_null_median_percentage_points_by_chain":[
                100*(obs_frac-a["median_fraction"]),
                100*(obs_frac-b["median_fraction"])],
            "null_evaluable_added_cell_mean_by_chain":[
                a["mean_null_added_cells"],b["mean_null_added_cells"]],
            "null_evaluable_added_cell_min_by_chain":[
                a["min_null_added_cells"],b["min_null_added_cells"]],
            "null_evaluable_added_cell_max_by_chain":[
                a["max_null_added_cells"],b["max_null_added_cells"]],
            "ESS_approx_by_chain":[a["ess_diagnostic"]["ess_approx"],
                                  b["ess_diagnostic"]["ess_approx"]],
            "nominal_serial_mcmc_right_tail_NOT_EXACT_P_by_chain":
                right_tails if ok else None,
            "qualified":ok,
        }
    return {
        "mode":mode,
        "status":"QUALIFIED_EXPLORATORY_PLANT_NETWORK_MCMC_NULL"
        if all(x["qualified"] for x in quality.values()) else
        "NOT_EVALUABLE_MCMC_OR_SWAP_SUPPORT",
        "qualification":quality,
        "comparisons":contrasts,
        "chain_receipts":chains,
    }


def run(source_root):
    model,regions,climate,raw,source=source_model(source_root)
    assert len(raw)==239 and len(regions)==365
    orig={h:set(cols) for h,cols in model["introduced_sets"].items()}
    assert sum(len(cols) for cols in orig.values())==EXPECTED_PLANT_EDGES
    masks,included=source_analog_masks(model,regions,climate,raw)
    observed,added_all=butterfly_metrics(model,orig,masks)
    assert added_all==14553 and len(included)==234
    for key,expected in EXPECTED_OBS.items():
        assert {k:observed[key][k] for k in expected}==expected,\
            f"Observed source mismatch: {key}: {observed[key]}"
    modes={}
    for mode in MODES:
        chains=[run_chain(model,regions,climate,orig,masks,mode,seed) for seed in SEEDS]
        modes[mode]=summarize_mode(chains,observed,mode)
    return {
        "schema":"chocho_20261019_host_network_preserving_climate_analogue_mcmc_v0.1",
        "status":"QUALIFIED_EXPLORATORY_ALL_MODES" if all(
            x["status"]=="QUALIFIED_EXPLORATORY_PLANT_NETWORK_MCMC_NULL"
            for x in modes.values()) else "PARTIAL_OR_NOT_EVALUABLE",
        "source_manifest_sha256":source["source_frozen_manifest_sha256"],
        "source_used":"Pinned original WCVP/HOSTS/CHELSA accepted-edge author-only archive",
        "species":239,"climate_evaluable_species":len(included),
        "plant_introduced_region_edges":EXPECTED_PLANT_EDGES,
        "original_butterfly_region_added_all":added_all,
        "original_butterfly_region_added_evaluable":14481,
        "draws_per_chain":DRAW_COUNT,"chain_seeds":list(SEEDS),
        "burn_attempts_per_plant_edge":BURN_FACTOR,
        "spacing_attempts_per_plant_edge":SPACING_FACTOR,
        "original_observed":observed,
        "modes":modes,
        "interpretation_limits":[
            "Original accepted butterfly-host identities and native plant ranges preserved; plant edges redistributed only within plant family and continent, with second mode additionally preserving coarse plant climatic distributions.",
            "Null DOES NOT preserve butterfly per-region added degrees or individual species butterfly gained area: only fraction of added opportunities is a comparable response across draws.",
            "Potential larval-resource geography, not documented local host use, adult butterfly survival, colonisation, population homogenisation or plant introduction history.",
            "MCMC serial-draw tails are NOT exact independent randomized p-values, and fit under source catalogues that can be incomplete.",
            "Coarse 2×2 climate quantiles do not preserve full plant climatic tolerance, biotic interactions or geographic autocorrelation.",
            "No raw licensed WCVP/HOSTS/TDWG dataset redistribution; derived aggregate estimates only."
        ]
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-root",type=Path,required=True)
    parser.add_argument("--output-json",type=Path,required=True)
    args=parser.parse_args()
    receipt=run(args.source_root)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",
                                encoding="utf-8")
    print(json.dumps({
        "status":receipt["status"],
        "original_observed":receipt["original_observed"],
        "modes":{
            key:{
                "status":val["status"],
                "qualification":val["qualification"],
                "comparisons":val["comparisons"],
                "burn_moved":[x["fraction_plant_edge_positions_changed_after_burn"]
                              for x in val["chain_receipts"]],
                "spacing_accepted":[x["accepted_swaps_per_spacing_mean"]
                              for x in val["chain_receipts"]],
            }
            for key,val in receipt["modes"].items()
        }
    },sort_keys=True))

#!/usr/bin/env python3
"""Two-chain plant-network-preserving posthoc climate-analogue resource null.

No invented butterfly outcomes: move observed introduced HOST PLANT×WGSRPD3
edges within plant Family×Level1×BIO1/BIO12 region climate bins. The fixed
butterfly-to-host graph and native plant ranges remain unchanged.
Changes to butterfly added-region totals are allowed; compare O/A, not O alone.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_butterfly_climate_analog_macro_20261019 import (
    source_model, qbounds
)
from analyze_climate_analog_fixedmargin_null_20261019 import ess_approx

QS=(.05,.10)
SEEDS=(20261019,20261020)
BURN_FACTOR=10
SPACING_FACTOR=4
N_DRAWS=199
EXPECTED={(0.05):(9451,14481), (0.10):(6948,14481)}


def signatures(model,regions,climate):
    covered=[i for i,v in enumerate(regions) if v in climate]
    masks={q:{} for q in QS}
    evaluated=[]
    for i,bits in enumerate(model["native_species"]):
        native=[j for j in covered if (bits>>j)&1]
        if len(native)<10:
            continue
        evaluated.append(i)
        for q in QS:
            lo1,hi1=qbounds([climate[regions[j]][0] for j in native],q)
            lo2,hi2=qbounds([climate[regions[j]][2] for j in native],q)
            mask=0
            for j in covered:
                record=climate[regions[j]]
                if lo1<=record[0]<=hi1 and lo2<=record[2]<=hi2:
                    mask|=1<<j
            masks[q][i]=mask
    assert len(evaluated)==234
    return evaluated,masks


def climate_eval(model,introduced,evaluated,masks):
    plant_bits={host:sum(1<<j for j in cells) for host,cells in introduced.items()}
    added_bits={}
    for i in evaluated:
        union=0
        for host in model["species_hosts"][i]:
            union |= plant_bits.get(host,0)
        added_bits[i]=union&~model["native_species"][i]
    A=sum(mask.bit_count() for mask in added_bits.values())
    if not A:raise ValueError("zero added resource denominator")
    return {
        str(q):{
            "O":sum((mask&masks[q][i]).bit_count()
                    for i,mask in added_bits.items()),
            "A":A,
            "share":sum((mask&masks[q][i]).bit_count()
                    for i,mask in added_bits.items())/A,
        }
        for q in QS
    }


def source_definitions(source_root,model,regions,climate):
    sys.path.insert(0,str(source_root/"code"))
    from analyze_butterfly_climate_matched_plant_range_null import (
        climate_strata, grouped_introduced_edges
    )
    classes,metadata=climate_strata(
        climate,regions,model["level1"],"temperature_precipitation_2x2")
    return classes,metadata,grouped_introduced_edges


def source_constraint_invariants(model,introduced):
    """Checksum-like semantic checks; no plant ID may gain/lose edges in a group."""
    total=sum(len(cells) for cells in introduced.values())
    assert total==30373, f"Unexpected plant×region edge budget: {total}"
    for host,cells in introduced.items():
        assert not (set(cells)&set(model["native_sets"].get(host,set())))
    return total


def run_chain(model,classes,evaluated,masks,grouping,seed):
    from analyze_butterfly_shared_host_structural_null import attempt_range_swaps
    intro={h:set(cells) for h,cells in model["introduced_sets"].items()}
    E=source_constraint_invariants(model,intro)
    orig={h:frozenset(cells) for h,cells in intro.items()}
    orig_per_plant={h:len(cells) for h,cells in intro.items()}
    groups,cumulative,diag=grouping(intro,model["host_families"],
                                    model["level1"],classes)
    assert diag["total_introduced_plant_region_edges"]==E
    if diag["climate_covered_fraction"]<.85 or diag["viable_fraction_of_total"]<.40:
        return {"status":"NOT_EVALUABLE_CLIMATE_SUPPORT_OR_SWAP_GEOMETRY",
                "seed":seed,"groups":diag}

    rng=random.Random(seed)
    accepted_burn=attempt_range_swaps(
        groups,cumulative,model["native_sets"],intro,BURN_FACTOR*E,rng)
    moved=sum(len(orig[h]^intro[h]) for h in intro)/(2*E)
    source_constraint_invariants(model,intro)
    assert all(len(intro[h])==n for h,n in orig_per_plant.items())
    counts={str(q):[] for q in QS}
    accepted=[]
    for i in range(N_DRAWS):
        acc=attempt_range_swaps(
            groups,cumulative,model["native_sets"],intro,SPACING_FACTOR*E,rng)
        accepted.append(acc)
        x=climate_eval(model,intro,evaluated,masks)
        for q in QS:
            counts[str(q)].append(x[str(q)])
        if i in (0,N_DRAWS-1):
            source_constraint_invariants(model,intro)
            assert all(len(intro[h])==n for h,n in orig_per_plant.items())
    metrics={}
    for q in QS:
        name=str(q)
        shares=[r["share"] for r in counts[name]]
        area=[r["A"] for r in counts[name]]
        analog=[r["O"] for r in counts[name]]
        metrics[name]={
            "share_stats":ess_approx(shares),
            "share_median":float(statistics.median(shares)),
            "added_resource_cells_median":float(statistics.median(area)),
            "added_resource_cells_range":[min(area),max(area)],
            "analog_cells_median":float(statistics.median(analog)),
            "null_share_history":shares,
            "added_region_total_history":area,
            "introduced_analog_count_history":analog,
        }
    return {
        "status":"DRAW_CHAIN_COMPLETE_UNQUALIFIED_UNTIL_PAIRED",
        "seed":seed,"E_introduced_plant_region_edges":E,
        "strata":diag,"burn_proposals":BURN_FACTOR*E,
        "burn_accepted":accepted_burn,
        "moved_fraction_from_original_after_burn":moved,
        "spacing_proposals":SPACING_FACTOR*E,
        "accepted_per_spacing_mean":float(statistics.mean(accepted)),
        "accepted_per_spacing_min":min(accepted),
        "n_draws":N_DRAWS,"by_q":metrics
    }


def run(source_root):
    model,regions,climate,rows,source=source_model(source_root)
    evaluated,masks=signatures(model,regions,climate)
    classes,classmeta,grouping=source_definitions(source_root,model,regions,climate)
    observed=climate_eval(model,model["introduced_sets"],evaluated,masks)
    assert all((observed[str(q)]["O"],observed[str(q)]["A"])==EXPECTED[q]
               for q in QS),observed
    original_edges=source_constraint_invariants(model,model["introduced_sets"])
    chain=[run_chain(model,classes,evaluated,masks,grouping,seed) for seed in SEEDS]
    if any(x["status"]!="DRAW_CHAIN_COMPLETE_UNQUALIFIED_UNTIL_PAIRED" for x in chain):
        return {"status":"NOT_EVALUABLE_CLIMATE_OR_GRAPH_SUPPORT",
                "observed":observed,"chains":chain}
    quality={};results={}
    for q in QS:
        key=str(q)
        stats=[x["by_q"][key]["share_stats"] for x in chain]
        median=[x["by_q"][key]["share_median"] for x in chain]
        ess=[v["ess_approx"] for v in stats]
        agree=abs(median[0]-median[1])<=2*max(v["std"] for v in stats)
        enough_ess=all(n>=30 for n in ess)
        swapped=all(x["moved_fraction_from_original_after_burn"]>=.10
                    for x in chain)
        accepted=all(x["burn_accepted"]>=100 and
                     x["accepted_per_spacing_mean"]>=100 for x in chain)
        good=agree and enough_ess and swapped and accepted
        quality[key]={
            "qualified":good,
            "between_chain_medians_agree":agree,
            "ess_both_ge_30":enough_ess,
            "burn_movement_both_ge_10pct":swapped,
            "accepted_swaps_ge_100":accepted,
        }
        actual=observed[key]["share"]
        tail=[(1+sum(a>=actual for a in c["by_q"][key]["null_share_history"]))/
              (len(c["by_q"][key]["null_share_history"])+1) for c in chain]
        results[key]={
            "observed":observed[key],
            "null_share_median_by_chain":median,
            "null_median_added_resources_by_chain":[
                c["by_q"][key]["added_resource_cells_median"] for c in chain],
            "observed_minus_null_share_pp_by_chain":[
                100*(actual-m) for m in median],
            "ess_by_chain":ess,
            "serial_right_tail_fractions_if_qualified":tail if good else None,
            "numerically_qualified":good,
            "null_changes_total_added_resource_cells":True,
        }

    return {
        "schema":"chocho_host_plant_family_continent_climate_fixed_network_null_v0.1",
        "status":(
            "QUALIFIED_EXPLORATORY_PLANT_NETWORK_CONDITIONED_NULL"
            if all(x["qualified"] for x in quality.values())
            else "NOT_EVALUABLE_PLANT_NULL_CHAIN_MIXING"
        ),
        "original_focal_butterflies":239,
        "eligible_native_host_climate_butterflies":len(evaluated),
        "introduced_plant_region_edges":original_edges,
        "climate_mode":"temperature_precipitation_2x2",
        "climate_strata_metadata":classmeta,
        "source_integrity":source,
        "butterfly_host_edges_kept_fixed":True,
        "per_plant_per_continent_coarse_climate_range_counts_fixed":True,
        "per_region_per_plant_family_introduced_counts_fixed":True,
        "native_host_plants_and_ranges_kept_fixed":True,
        "newly_added_butterfly_region_count_changes_under_null":True,
        "seeds":SEEDS,
        "burn_factor":BURN_FACTOR,
        "spacing_factor":SPACING_FACTOR,
        "draws_per_chain":N_DRAWS,
        "quality":quality,
        "result":results,
        "chains":chain,
        "limits":[
            "Posthoc model chosen after observing positive geometry-conditioned climate alignment.",
            "Coarse BIO1/BIO12 climate bins leave within-bin climatic variation available for swap.",
            "This is climate-analog potential plant resource geography, NOT butterfly occupancy, reproduction or dispersal.",
            "The 199 MCMC draws are serial dependent; add-one tails are not exact permutation p-values.",
            "Effects from the incomplete HOSTS butterfly-host catalogue cannot establish an ecosystem-wide evolutionary rule.",
            "No current source licenses to redistribute author-only WCVP/TDWG source archive."
        ],
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    ans=run(args.source_root)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(ans,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":ans["status"],
        "quality":ans.get("quality"),
        "result":ans.get("result"),
        "moved_after_burn":[c.get("moved_fraction_from_original_after_burn") for c in ans["chains"]],
        "viable_edges":[c.get("strata",{}).get("viable_edges") for c in ans["chains"]],
    },sort_keys=True))

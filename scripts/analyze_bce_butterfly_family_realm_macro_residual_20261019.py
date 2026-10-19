#!/usr/bin/env python3
"""239 butterfly resource-climate residual partition by taxonomic family and Level1.

Source-verified, grouped PRIMARY ecological-heterogeneity analysis of the
already-qualified accepted plant-family × Level1 × BIO1/BIO12 3x3 null.
Do not infer butterfly colonisation or any insect survival/recruitment.
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
    source_augmented_model, analog_mask, marginal_signature,
    EXPECTED_INTRODUCED_PLANT_REGION_EDGES, SEEDS,
)
from analyze_bce_finer_climate_spatial_null_20261019 import source_stratum_setup
from analyze_climate_analog_fixedmargin_null_20261019 import ess_approx
from analyze_global_butterfly_habitat_filter_20261019 import read_inputs

MODE = "plant_bio1_bio12_3x3"
COHORTS = ("augmented_235", "matched_original_234")
Q = (0.05, 0.10)
SOURCES = {
    "augmented_235": {"0.05":(9949,14961),"0.1":(7370,14961)},
    "matched_original_234":{"0.05":(9903,14903),"0.1":(7331,14903)},
}
SEED_N = 199
BURN_FACTOR = 10
SPACING_FACTOR = 4
MIN_SPECIES_FAMILY = 5
MIN_ADDED = 150
MIN_ESS = 30


def family_index(model):
    rows, _, hashes = read_inputs()
    species_family = {r["species"]:r["family"] for r in rows}
    assert len(species_family)==239
    if set(model["species"]) != set(species_family):
        raise ValueError("Butterfly family source and accepted host species disagree")
    family_for_i = {i:species_family[species] for i,species in enumerate(model["species"])}
    assert len(set(family_for_i.values()))==6
    return family_for_i, hashes


def realm_index(model):
    groups=defaultdict(int)
    for idx,region in model["level1"].items():
        groups[str(region)] |= (1<<idx)
    assert len(groups)==9
    return dict(sorted(groups.items()))


def grouped_opportunity(model, plant_sets, analogs, fams, realms):
    # A is never held fixed in the plant swap counterfactual.
    host_bits = {h:sum(1 << j for j in cols) for h,cols in plant_sets.items()}
    added = []
    for i,links in enumerate(model["species_hosts"]):
        union=0
        for h in links:
            union |= host_bits.get(h,0)
        added.append(union & ~model["native_species"][i])
    results={}
    for cohort in COHORTS:
        for q in Q:
            eligible=analogs[q]
            eligible_indices=[
                i for i in eligible
                if cohort=="augmented_235" or model["species"][i]!="Satyrium pruni"
            ]
            all_groups={
                "all":{"all":{"O":0,"A":0,"n_species":0,"n_added_species":0}},
                "family":{k:{"O":0,"A":0,"n_species":0,"n_added_species":0}
                           for k in sorted(set(fams.values()))},
                "realm":{k:{"O":0,"A":0,"n_species":0,"n_added_species":0}
                           for k in realms},
            }
            # Precompute integer mask intersections. Region receivers can overlap
            # taxonomically but are disjoint geographical Level1 strata.
            for i in eligible_indices:
                mask=added[i]
                aligned=mask & eligible[i]
                N=mask.bit_count()
                O=aligned.bit_count()
                f=fams[i]
                for target in (all_groups["all"]["all"],all_groups["family"][f]):
                    target["O"]+=O
                    target["A"]+=N
                    target["n_species"]+=1
                    target["n_added_species"]+=int(N>0)
                for k,bits in realms.items():
                    a=(mask & bits).bit_count()
                    o=(aligned & bits).bit_count()
                    target=all_groups["realm"][k]
                    target["O"]+=o
                    target["A"]+=a
                    if a:
                        target["n_added_species"]+=1
            for level,table in all_groups.items():
                for g,datum in table.items():
                    datum["fraction"]=datum["O"]/datum["A"] if datum["A"] else None
            results[cohort,str(q)]=all_groups
    return results


def check_source_observed(baseline):
    for cohort,tests in SOURCES.items():
        for q,(o,a) in tests.items():
            observed=baseline[cohort,q]["all"]["all"]
            if (observed["O"],observed["A"]) != (o,a):
                raise ValueError(f"Source-baseline drift {cohort} q={q}: {observed}")
    for cohort in COHORTS:
        for q in map(str,Q):
            r=baseline[cohort,q]
            A=r["all"]["all"]["A"];O=r["all"]["all"]["O"]
            assert sum(x["A"] for x in r["family"].values()) == A
            assert sum(x["O"] for x in r["family"].values()) == O
            assert sum(x["A"] for x in r["realm"].values()) == A
            assert sum(x["O"] for x in r["realm"].values()) == O
    assert baseline["augmented_235","0.05"]["all"]["all"]["n_species"]==235
    assert baseline["matched_original_234","0.05"]["all"]["all"]["n_species"]==234


def run_chain(model,regions,climate,analog,original,fams,realms,seed,
              *, n=SEED_N, burn=BURN_FACTOR, spacing=SPACING_FACTOR):
    sys.path.insert(0,str(SOURCE_ROOT/"code"))
    from analyze_butterfly_shared_host_structural_null import attempt_range_swaps
    plant_sets={h:set(cells) for h,cells in original.items()}
    classes,groups,cumulative,structure,metadata,spatial=source_stratum_setup(
        SOURCE_ROOT,model,regions,climate,plant_sets,MODE)
    if spatial:
        raise ValueError("Unexpected spatial group in prespecified 3x3 model")
    total=sum(map(len,plant_sets.values()))
    assert total==EXPECTED_INTRODUCED_PLANT_REGION_EDGES==37869
    assert structure["viable_fraction_of_total"]>=.75
    baseline=marginal_signature(model,original,classes)
    assert marginal_signature(model,plant_sets,classes)==baseline
    rng=random.Random(seed)
    accepted_burn=attempt_range_swaps(
        groups,cumulative,model["native_sets"],plant_sets,burn*total,rng)
    changed=sum(len(plant_sets[h]^original[h]) for h in original)/(2*total)
    assert marginal_signature(model,plant_sets,classes)==baseline

    data={
        cohort:{str(q):{g:{} for g in ("all","family","realm")} for q in Q}
        for cohort in COHORTS
    }
    accepted=[]
    for t in range(n):
        accepted.append(attempt_range_swaps(
            groups,cumulative,model["native_sets"],plant_sets,spacing*total,rng))
        # One plant graph -> all groups, both cohorts and both native-host
        # climate analogue widths. No independent group-wise simulations.
        group=grouped_opportunity(model,plant_sets,analog,fams,realms)
        for cohort,q in group:
            for level,table in group[cohort,q].items():
                for label,row in table.items():
                    history=data[cohort][q][level].setdefault(
                        label,{"fractions":[],"denominators":[]})
                    history["fractions"].append(row["fraction"])
                    history["denominators"].append(row["A"])
        if t in (0,n//2,n-1):
            if marginal_signature(model,plant_sets,classes)!=baseline:
                raise RuntimeError("Family, stratum, native or region degree drift")

    output={}
    for cohort in COHORTS:
        output[cohort]={}
        for q in map(str,Q):
            output[cohort][q]={}
            for level,table in data[cohort][q].items():
                output[cohort][q][level]={}
                for label,record in table.items():
                    if len(record["fractions"])!=n:raise RuntimeError("incomplete group draw")
                    complete=all(v is not None for v in record["fractions"])
                    if complete:
                        x=np.array(record["fractions"],dtype=float)
                        ess=ess_approx(x.tolist())
                        med=float(np.median(x))
                        ci=[float(np.quantile(x,z)) for z in (.025,.975)]
                        null_std=float(np.std(x,ddof=1))
                    else:
                        ess={"ess_approx":0,"reason":"Some draw denominators zero"}
                        med=None;ci=None;null_std=None
                    output[cohort][q][level][label]={
                        "null_median":med,
                        "null_percentile_95":ci,
                        "null_std":null_std,
                        "ess_diagnostic":ess,
                        "null_denominator_median":float(np.median(record["denominators"])),
                        "null_denominator_min":int(min(record["denominators"])),
                        "null_denominator_max":int(max(record["denominators"])),
                        "draw_count":n,
                    }
    return {
        "seed":seed,"mode":MODE,"sampled":n,
        "accepted_burn":accepted_burn,
        "moved_fraction":changed,
        "mean_accepted_per_spacing":float(np.mean(accepted)),
        "source_viable_edge_fraction":structure["viable_fraction_of_total"],
        "original_network_margin_checks":"After burn + first/middle/last sampling draws",
        "summary":output,
    }


def merge_group_results(observed,chain,*,global_gate):
    result={}
    for cohort in COHORTS:
        result[cohort]={}
        for q in map(str,Q):
            result[cohort][q]={}
            for level,groups in observed[cohort,q].items():
                result[cohort][q][level]={}
                for label,raw in groups.items():
                    n=raw["n_added_species"]
                    A=raw["A"]
                    support=(A>=MIN_ADDED and
                             (level!="family" or n>=MIN_SPECIES_FAMILY))\
                            if level!="all" else True
                    ss=[c["summary"][cohort][q][level][label] for c in chain]
                    nullmed=[x["null_median"] for x in ss]
                    ess=[x["ess_diagnostic"].get("ess_approx",0) for x in ss]
                    valid_med=all(x is not None for x in nullmed)
                    sds=[x["null_std"] for x in ss]
                    med_agree=(valid_med and all(x is not None for x in sds)
                               and abs(nullmed[0]-nullmed[1])<=2*max(sds))
                    quality=(global_gate and support and valid_med and med_agree
                             and all(e>=MIN_ESS for e in ess))
                    result[cohort][q][level][label]={
                        "status":("QUALIFIED_EXPLORATORY_HETEROGENEITY" if quality
                                  else "NOT_EVALUABLE_GROUP_SUPPORT_OR_MCMC"),
                        "species_with_added":n,
                        "species_total":raw["n_species"] if level!="realm" else None,
                        "observed_aligned":raw["O"],
                        "observed_added_cells":A,
                        "observed_alignment_fraction":raw["fraction"],
                        "null_fraction_median_by_chain":nullmed,
                        "null_fraction_95_by_chain":[x["null_percentile_95"] for x in ss],
                        "null_added_median_by_chain":[x["null_denominator_median"] for x in ss],
                        "ESS_approx_by_chain":ess,
                        "observed_minus_null_median_pp_by_chain":[
                            100*(raw["fraction"]-med) for med in nullmed
                        ] if quality else None,
                        "adequate_group_support":bool(support),
                        "median_agreement":bool(med_agree),
                        "full_global_MCMC_structure_pass":bool(global_gate),
                    }
    return result


SOURCE_ROOT=Path("unused")

def run(source_root,*,n=SEED_N):
    global SOURCE_ROOT
    SOURCE_ROOT=Path(source_root)
    model,regions,climate,source=source_augmented_model(SOURCE_ROOT)
    masks,eligible=analog_mask(model,regions,climate)
    assert len(eligible)==235 and len(regions)==365
    fams,input_hashes=family_index(model)
    realms=realm_index(model)
    original={h:set(cells) for h,cells in model["introduced_sets"].items()}
    observed=grouped_opportunity(model,original,masks,fams,realms)
    check_source_observed(observed)
    chain=[run_chain(model,regions,climate,masks,original,fams,realms,seed,n=n)
           for seed in SEEDS]
    global_gate=all(c["moved_fraction"]>=.10 and
                    c["mean_accepted_per_spacing"]>=100 and
                    c["source_viable_edge_fraction"]>=.75 for c in chain)
    summaries=merge_group_results(observed,chain,global_gate=global_gate)
    global_ess_pass=all(
        summaries[cohort][q]["all"]["all"]["status"]==
        "QUALIFIED_EXPLORATORY_HETEROGENEITY"
        for cohort in COHORTS for q in map(str,Q))
    if not global_ess_pass:
        status="NOT_EVALUABLE_GLOBAL_MCMC"
    else:
        status="GLOBAL_SOURCE_QUALIFIED_GROUP_HETEROGENEITY_EXPLORED"
    return {
        "schema":"chocho_20261019_239butterfly_family_and_level1_climate_residual_v0.1",
        "status":status,
        "source_provenance":source,
        "family_data_source_digests":input_hashes,
        "mcmc_mode":MODE,
        "sources":{
            "239_focal_species":len(model["species"]),
            "butterfly_families":sorted(set(fams.values())),
            "recipient_level1_codes":sorted(realms),
            "total_introduced_plant_region_edges":sum(map(len,original.values())),
            "source_augmented_eligible":235,"source_matched_eligible":234,
            "reference_wide_added":14961,"reference_narrow_added":14961
        },
        "chain_quality":[{
            "seed":c["seed"],"accepted_burn":c["accepted_burn"],
            "moved_fraction":c["moved_fraction"],
            "mean_accepted_per_spacing":c["mean_accepted_per_spacing"],
            "source_viable_edge_fraction":c["source_viable_edge_fraction"],
        } for c in chain],
        "family_and_realm_results":summaries,
        "bounds":[
            "Butterfly taxonomic-family and WGSRPD Level1 grouping are not independent ecological experiments; butterflies may share hosts and source areas.",
            "Species identification and source host database incompleteness are major biases; BCE catalogue additions are Europe-biased.",
            "Small or MCMC-nonconverged groups cannot support directional ecological generalization.",
            "All outcomes are potential larval host resource climate analog fractions; not observed butterfly occupation, life-cycle survival, oviposition or population change.",
            "WGSRPD3 cells have unequal physical areas, and original WCVP/HOSTS provider raw geography is not redistributed.",
        ]
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--draws",type=int,default=SEED_N)
    args=ap.parse_args()
    ans=run(args.source_root,n=args.draws)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(ans,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":ans["status"],
        "source_239":ans["sources"],
        "family_results":{
            k:{family:{
                "status":v["status"],
                "observed_added":v["observed_added_cells"],
                "observed_fraction":v["observed_alignment_fraction"],
                "residual_pp":v["observed_minus_null_median_pp_by_chain"],
                "ESS":v["ESS_approx_by_chain"]
            } for family,v in ans["family_and_realm_results"]["augmented_235"][k]["family"].items()}
            for k in ("0.05","0.1")
        },
        "realm_results":{
            k:{realm:{
                "status":v["status"],
                "observed_added":v["observed_added_cells"],
                "residual_pp":v["observed_minus_null_median_pp_by_chain"],
                "ESS":v["ESS_approx_by_chain"]
            } for realm,v in ans["family_and_realm_results"]["augmented_235"][k]["realm"].items()}
            for k in ("0.05","0.1")
        }
    },sort_keys=True))

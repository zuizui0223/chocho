#!/usr/bin/env python3
"""Macroecological climate-analogue test using pinned HOSTS/WCVP/CHELSA source.

Original source archive is retrieved locally by Github Actions from a pre-existing
author-only CI artifact; extracted contents are NOT redistributed.
No butterfly occupancy, colonisation, phenology or larval survival is measured.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
from analyze_global_butterfly_habitat_filter_20261019 import (
    read_inputs, annotate, adjusted, STRONG_ASSOC, STRONG_AVOID
)

EXPECTED_INNER_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
SEED=20261019

def read_csv(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def source_model(source:Path):
    manifest=json.loads((source/"SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["schema"]=="butterfly_anonymous_offline_resource_package_v0.1"
    assert manifest["inputs_count"]==11
    verified=0
    for name,digest in manifest["sha256"].items():
        path=source/name
        assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==digest, name
        verified+=1
    sys.path.insert(0,str(source/"code"))
    from analyze_butterfly_shared_host_structural_null import build_model
    from analyze_butterfly_climate_matched_plant_range_null import read_unit_climate
    data=source/"inputs"
    all_hosts=read_csv(data/"original/insect_host_accepted.csv")
    native=read_csv(data/"original/native_extant_nondoubtful_wgsrpd3.csv")
    contemporary=read_csv(data/"original/extant_nondoubtful_wgsrpd3.csv")
    desc=read_csv(data/"s1_resource_descriptors.csv")
    lookup=json.loads((data/"level1_lookup.json").read_text(encoding="utf-8"))["units"]
    model=build_model(desc,all_hosts,native,contemporary,lookup)
    assert len(model["species"])==239
    all_used_hosts={host for hs in model["species_hosts"] for host in hs}
    regions=sorted({
        row["area_code_l3"] for row in native+contemporary
        if row["accepted_plant_name_id"] in all_used_hosts
    })
    assert len(regions)==365 and set(regions)<=set(lookup)
    clim,provenance=read_unit_climate(data/"unit_climate_crossfit.csv")
    assert provenance["source_unique_valid_regions"]>=300
    climate_codes=set(clim)
    covered_regions={i for i,r in enumerate(regions) if r in climate_codes}
    assert len(covered_regions)>=300
    source_metric=read_csv(data/"anthropogenic_species_metrics.csv")
    metric={r["species"]:r for r in source_metric}
    assert len(metric)==239
    ranges={h:sum(1 << r for r in cells) for h,cells in model["introduced_sets"].items()}
    data_rows=[]
    for i,s in enumerate(model["species"]):
        native_bits=model["native_species"][i]
        introduced_bits=0
        for host in model["species_hosts"][i]:
            introduced_bits|=ranges.get(host,0)
        added_bits=introduced_bits & ~native_bits
        n=native_bits.bit_count()
        a=added_bits.bit_count()
        c=(native_bits|added_bits).bit_count()
        original=metric[s]
        assert n==int(original["native_resource_units"]), ("native",s,n)
        assert c==int(original["contemporary_resource_units"]),("current",s,c)
        assert a==int(original["introduced_added_units"]),("added",s,a)
        data_rows.append({
            "species":s,"genus":s.split(" ")[0],
            "native_bits":native_bits,"added_bits":added_bits,
            "native":n,"added":a,"contemporary":c
        })
    assert sum(r["native"] for r in data_rows)==26530
    assert sum(r["added"] for r in data_rows)==14553
    return model,regions,clim,data_rows,{
        "source_members_verified":verified,
        "source_frozen_manifest_sha256":hashlib.sha256(
            (source/"SOURCE_MANIFEST.json").read_bytes()).hexdigest(),
        "climate":provenance,
        "mode":"AUTHOR_ONLY_RESTORED_SOURCE_NO_REDISTRIBUTION"
    }


def bits_to_indices(mask,n):
    return [i for i in range(n) if mask >> i & 1]


def qbounds(values,lower):
    vals=np.array(values,dtype=float)
    return float(np.quantile(vals,lower)),float(np.quantile(vals,1-lower))


def calculate_analog(model,regions,climate,raw_rows,q):
    N=len(regions)
    covered=[i for i,s in enumerate(regions) if s in climate]
    annual_temp={i:float(climate[regions[i]][0]) for i in covered}
    annual_precip={i:float(climate[regions[i]][2]) for i in covered}
    by_l1=defaultdict(list)
    for i in covered:
        by_l1[model["level1"][i]].append(i)
    result=[]
    skip=Counter()
    for raw in raw_rows:
        native={i for i in bits_to_indices(raw["native_bits"],N)}
        introduced={i for i in bits_to_indices(raw["added_bits"],N)}
        native_clim=native&set(covered)
        introduced_clim=introduced&set(covered)
        if len(native_clim)<10:
            skip["native_climate_support_lt_10"]+=1
            continue
        if len(introduced-introduced_clim)>0:
            skip["introduced_added_missing_any_climate"]+=1
            continue
        temp_lo,temp_hi=qbounds([annual_temp[i] for i in native_clim],q)
        precip_lo,precip_hi=qbounds([annual_precip[i] for i in native_clim],q)
        analog={
            i for i in covered if temp_lo<=annual_temp[i]<=temp_hi
             and precip_lo<=annual_precip[i]<=precip_hi
        }
        cap=set(covered)-native
        eligible_analog=cap&analog
        obs=int(len(introduced_clim&analog))
        expected=0.0
        reference_var=0.0
        added_by_l1=Counter(model["level1"][i] for i in introduced_clim)
        for group,indices in by_l1.items():
            candidate=set(indices)-native
            if not candidate:
                assert added_by_l1[group]==0
                continue
            M=len(candidate)
            k=len(candidate&analog)
            n_add=added_by_l1[group]
            assert 0<=n_add<=M
            p=k/M
            expected+=n_add*p
            if M>1:
                reference_var+=n_add*p*(1-p)*(M-n_add)/(M-1)
        K=len(eligible_analog)
        A=len(introduced_clim)
        assert obs<=K and obs<=A
        result.append({
            "species":raw["species"],"genus":raw["genus"],
            "native_regions":raw["native"],"added_regions":raw["added"],
            "native_climate_supported_regions":len(native_clim),
            "climate_covered_added_regions":A,
            "all_climate_covered_candidate_non_native_regions":len(cap),
            "eligible_climate_analog_non_native_regions":K,
            "observed_added_analog_regions":obs,
            "within_continent_reference_expected_analog_regions":expected,
            "within_continent_reference_variance":reference_var,
            "observed_analog_share_among_added":obs/A if A else None,
            "reference_analog_share_among_added":expected/A if A else None,
            "analog_excess_share":(obs-expected)/A if A else None,
            "analog_fill_fraction":obs/K if K else None,
            "analog_fill_logit":math.log((obs+0.5)/(K-obs+0.5)) if K else None,
            "climate_q_low":q
        })
    return result,dict(sorted(skip.items()))


def bootstrap_genus(rows,*,seed=SEED,n=2499):
    members=defaultdict(list)
    for r in rows:
        members[r["genus"]].append(r)
    genera=sorted(members)
    rng=np.random.default_rng(seed)
    samples=[]
    for i in range(n):
        sample=[r for j in rng.integers(0,len(genera),len(genera))
                for r in members[genera[int(j)]]]
        mean_delta=float(np.mean([r["analog_excess_share"] for r in sample]))
        obs=sum(r["observed_added_analog_regions"] for r in sample)
        added=sum(r["added_regions"] for r in sample)
        exp=sum(r["within_continent_reference_expected_analog_regions"] for r in sample)
        samples.append((mean_delta,(obs-exp)/added if added else None))
    a=np.array([x[0] for x in samples])
    b=np.array([x[1] for x in samples if x[1] is not None])
    return {
        "resampled_genera":len(genera),
        "draws":n,
        "species_equal_mean_excess_share_95ci":[float(np.quantile(a,q)) for q in (.025,.975)],
        "added_cell_weighted_excess_share_95ci":[float(np.quantile(b,q)) for q in (.025,.975)]
    }


def summarise_analog(rows,exclusions):
    usable=[r for r in rows if r["added_regions"]>0]
    tot_A=sum(r["added_regions"] for r in rows)
    obs=sum(r["observed_added_analog_regions"] for r in rows)
    exp=sum(r["within_continent_reference_expected_analog_regions"] for r in rows)
    delta=[r["analog_excess_share"] for r in usable]
    assert tot_A>0
    return {
        "species_with_adequate_native_climate":len(rows),
        "with_any_introduced_added":len(usable),
        "species_zero_added":len(rows)-len(usable),
        "species_with_zero_eligible_analog_candidates":sum(r["eligible_climate_analog_non_native_regions"]==0 for r in rows),
        "added_cells_all":tot_A,
        "observed_added_climate_analog":obs,
        "within_continent_reference_expected_analog":exp,
        "observed_analog_fraction_of_added":obs/tot_A,
        "expected_analog_fraction_of_added":exp/tot_A,
        "added_cell_weighted_excess_fraction":(obs-exp)/tot_A,
        "species_equal_excess_fraction":float(np.mean(delta)),
        "species_equal_fraction_with_observed_gt_expected":sum(x>1e-10 for x in delta)/len(delta),
        "species_equal_fraction_with_observed_lt_expected":sum(x<-1e-10 for x in delta)/len(delta),
        "genera_bootstrap":bootstrap_genus(usable),
        "missing_support":exclusions,
        "benchmark_limitation":"Within-WGSRPD Level1 uniform cell allocation fixes butterfly-added cells per broad area, not plant distributions or local climate suitability; NOT a causal/null-permutation p-value."
    }


def trait_climate_test(clim_rows,macro_rows, *, q, minK=5):
    mapped={r["species"]:r for r in clim_rows}
    filtered=[]
    for r in annotate(macro_rows):
        c=mapped.get(r["species"])
        if c is None or c["eligible_climate_analog_non_native_regions"]<minK:
            continue
        x=r.copy()
        x["climate_fill_logit"]=c["analog_fill_logit"]
        filtered.append(x)
    n_assoc=sum(x["exposure"] for x in filtered)
    n_avoid=len(filtered)-n_assoc
    out={
        "eligible_species":len(filtered),
        "n_disturbance_associated_strong":n_assoc,
        "n_disturbance_avoidant_strong":n_avoid,
        "min_eligible_climate_analogue_regions":minK,
        "climate_percentile_lower":q,
        "response":"continuity_correction_logit[(introduced_added_climate_analog+0.5)/(K-introduced_added_climate_analog+0.5)]",
        "zero_added_species_retained":sum(x["added"]==0 for x in filtered),
        "status":None
    }
    if len(filtered)<80 or n_avoid<12 or n_assoc<30:
        out["status"]="NOT_EVALUABLE_TRAIT_CLASS_SUPPORT"
        return out
    all_families=sorted({r["family"] for r in macro_rows})
    test=adjusted(filtered,all_families,
       label="strict_climate_analog_capacity_filter",
       y_key="climate_fill_logit",bootstrap=True)
    out.update({"status":"EXPLORATORY_ADJUSTED_ASSOCIATION","estimate":test})
    return out


def run(source_root):
    model,regions,climate,raw,source=source_model(source_root)
    macro,traits,digest=read_inputs()
    assert len(raw)==239 and len(macro)==239
    result={}
    for q in (.05,.10):
        rows,excluded=calculate_analog(model,regions,climate,raw,q)
        result[str(q)]={
            "overview":summarise_analog(rows,excluded),
            "secondary_disturbance_5_region_gate":trait_climate_test(rows,macro,q=q,minK=5),
            "secondary_disturbance_10_region_gate":trait_climate_test(rows,macro,q=q,minK=10),
            "species_derived":sorted(rows,key=lambda x:x["species"])
        }
    return {
        "schema":"chocho_global_butterfly_climate_analog_opportunity_v0.1",
        "status":"EXPLORATORY_POTENTIAL_HOST_RESOURCE_CLIMATE_FILTER_NOT_REALIZED_BUTTERFLY_OCCUPANCY",
        "scope":"239 focal butterflies across 365 WGSRPD3 resource regions",
        "source_integrity":source,
        "macro_input_digests":digest,
        "climate_variables":["CHELSA bio1","CHELSA bio12"],
        "main_threshold_5_to_95":result["0.05"],
        "sensitivity_threshold_10_to_90":result["0.1"],
        "negative_gates":{
            "adult_butterfly_recruitment_measured":False,
            "butterfly_native_climatic_physiology_measured":False,
            "native_host_climate_proxy_only":True,
            "within_continent_expectation_preserves_plant_spatial_network":False,
            "region_cells_have_equal_area":False,
            "source_archive_redistribution_authorised":False,
            "matched_butterfly_occupation_and_historical_host_plant_introduction":False,
            "confirmatory_ecological_causation":False,
        }
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    v=run(args.source_root)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    report={}
    for section in ("main_threshold_5_to_95","sensitivity_threshold_10_to_90"):
        d=v[section]
        report[section]={
            "overview":{k:z for k,z in d["overview"].items() if k not in ("benchmark_limitation",)},
            "disturbance_K5":d["secondary_disturbance_5_region_gate"]["status"],
            "disturbance_K5_n":d["secondary_disturbance_5_region_gate"]["eligible_species"],
            "disturbance_K5_assoc":d["secondary_disturbance_5_region_gate"]["n_disturbance_associated_strong"],
            "disturbance_K5_avoid":d["secondary_disturbance_5_region_gate"]["n_disturbance_avoidant_strong"],
            "disturbance_K5_beta":d["secondary_disturbance_5_region_gate"].get("estimate",{}).get("adjusted_coefficient"),
            "disturbance_K5_cluster_ci":d["secondary_disturbance_5_region_gate"].get("estimate",{}).get("approx_genus_cluster_95ci"),
            "disturbance_K5_bootstrap_ci":d["secondary_disturbance_5_region_gate"].get("estimate",{}).get("genus_bootstrap",{}).get("percentile_95ci"),
            "disturbance_K10_status":d["secondary_disturbance_10_region_gate"]["status"],
            "disturbance_K10_beta":d["secondary_disturbance_10_region_gate"].get("estimate",{}).get("adjusted_coefficient")
        }
    print(json.dumps(report,sort_keys=True))

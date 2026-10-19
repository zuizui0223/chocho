#!/usr/bin/env python3
"""Butterfly×WGSRPD3 fixed-margin climate alignment null, original 239 species.

This is a posthoc exploratory stricter structural counterfactual, NOT an exact
plant-network randomization or observed butterfly colonisation.
Requires the same author-only, source-byte-verified offline inputs as
analyze_butterfly_climate_analog_macro_20261019.py.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from analyze_butterfly_climate_analog_macro_20261019 import (
    source_model, bits_to_indices, qbounds
)

SEEDS=(20261019,20261020)
EXPECTED={"0.05":9451,"0.1":6948}
Q=(.05,.10)
BURN_FACTOR=10
SPACING_FACTOR=4
DRAW_COUNT=199


def analog_masks(model,regions,climate,raw):
    covered={i for i,code in enumerate(regions) if code in climate}
    out={q:{} for q in Q}
    included=[]
    for i,row in enumerate(raw):
        native=set(bits_to_indices(row["native_bits"],len(regions)))
        eligible=sorted(native&covered)
        if len(eligible)<10:
            continue
        included.append(i)
        for q in Q:
            lo1,hi1=qbounds([climate[regions[j]][0] for j in eligible],q)
            lo2,hi2=qbounds([climate[regions[j]][2] for j in eligible],q)
            mask=0
            for j in covered:
                a=climate[regions[j]]
                if lo1<=a[0]<=hi1 and lo2<=a[2]<=hi2:
                    mask |= 1<<j
            out[q][i]=mask
    assert len(included)==234
    return out,included


def count_alignment(species_added_bits, analog_masks):
    return {
        str(q):sum((species_added_bits[i]&m).bit_count() for i,m in analog_masks[q].items())
        for q in Q
    }


def create_area_graph(model,raw):
    groups=defaultdict(list)
    native=[x["native_bits"] for x in raw]
    added_original=[x["added_bits"] for x in raw]
    for b,row in enumerate(raw):
        for region in bits_to_indices(row["added_bits"],365):
            groups[model["level1"][region]].append([b,region])
    groups={k:v for k,v in sorted(groups.items())}
    E=sum(len(v) for v in groups.values())
    assert E==14553 and sum(m.bit_count() for m in added_original)==E
    return groups,native,added_original


def check_margins(groups,original_native,orig_added,baseline_region_degree):
    added=[0]*len(orig_added)
    regional=Counter()
    per_butterfly_area=Counter()
    for key,edges in groups.items():
        for b,region in edges:
            assert not ((original_native[b]>>region)&1)
            assert not ((added[b]>>region)&1)
            added[b] |= 1<<region
            regional[region]+=1
            per_butterfly_area[b,key]+=1
    assert sum(x.bit_count() for x in added)==14553
    assert regional==baseline_region_degree
    original_per_butterfly_area=Counter()
    for b,mask in enumerate(orig_added):
        for i in bits_to_indices(mask,365):
            original_per_butterfly_area[b,_LEVEL1_BY_REGION[i]]+=1
    assert per_butterfly_area==original_per_butterfly_area
    return added


# Bound during each run; never global cached across workflows.
_LEVEL1_BY_REGION={}


def attempt(groups,keys,counts,cumulative,added_bits,native,rng,n_attempts):
    accepted=0
    for _ in range(n_attempts):
        # Select a focal added-incidence at uniform probability, then another
        # within its already-fixed Level1 stratum. Proposal is symmetric.
        k=rng.randrange(cumulative[-1])
        group_idx=bisect.bisect_right(cumulative,k)
        lo=cumulative[group_idx-1] if group_idx else 0
        edges=groups[keys[group_idx]]
        i=k-lo
        j=rng.randrange(len(edges))
        if i==j:
            continue
        b1,r1=edges[i]
        b2,r2=edges[j]
        if b1==b2 or r1==r2:
            continue
        if ((native[b1]>>r2)&1) or ((native[b2]>>r1)&1):
            continue
        if ((added_bits[b1]>>r2)&1) or ((added_bits[b2]>>r1)&1):
            continue
        # Swap ONLY butterfly targets; region counts & butterfly-by-Level1
        # margins stay fixed and existing native-host cells remain unavailable.
        added_bits[b1] &= ~(1<<r1)
        added_bits[b1] |= 1<<r2
        added_bits[b2] &= ~(1<<r2)
        added_bits[b2] |= 1<<r1
        edges[i][1]=r2
        edges[j][1]=r1
        accepted+=1
    return accepted


def ess_approx(draws):
    a=np.asarray(draws,dtype=float)
    n=len(a)
    if n<20 or np.var(a)<=1e-20:
        return {"n":n,"lag1":None,"ess_approx":0.0,"mean":float(np.mean(a)),
                "std":float(np.std(a)),"early_late_drift":None}
    x=a-np.mean(a)
    den=float(x@x)
    coeffs=[]
    for lag in range(1,min(n//3,30)):
        rho=float(x[:-lag]@x[lag:]/den) # biased autocorr denominator for stability
        if rho<=0: break
        coeffs.append(rho)
    iact=1+2*sum(coeffs)
    return {
        "n":n,"lag1":float(x[:-1]@x[1:]/den),
        "ess_approx":float(n/iact),
        "iact_positive_seq":iact,
        "mean":float(np.mean(a)),
        "median":float(np.median(a)),
        "std":float(np.std(a,ddof=1)),
        "early_late_drift":float(np.mean(a[-n//4:])-np.mean(a[:n//4])),
    }


def run_chain(model,raw,analog,seed):
    groups,native,initial=create_area_graph(model,raw)
    E=sum(len(x) for x in groups.values())
    added=initial[:]
    keys=sorted(groups)
    counts=[len(groups[k]) for k in keys]
    cumulative=[]
    s=0
    for c in counts:
        s+=c
        cumulative.append(s)
    assert s==E
    rng=random.Random(seed)
    burn=attempt(groups,keys,counts,cumulative,added,native,rng,BURN_FACTOR*E)
    moved=sum((initial[i]&~added[i]).bit_count() for i in range(239))/E
    originals=Counter(r for x in raw for r in bits_to_indices(x["added_bits"],365))
    samples={str(q):[] for q in Q}
    accepted=[]
    for k in range(DRAW_COUNT):
        n=attempt(groups,keys,counts,cumulative,added,native,rng,SPACING_FACTOR*E)
        accepted.append(n)
        if k in (0,DRAW_COUNT-1):
            # Rebuild native-excluded baseline and every region degree; test
            # species×continent counts separately from chain state.
            check_margins(groups,native,initial,originals)
        obs=count_alignment(added,analog)
        for q in Q:
            samples[str(q)].append(obs[str(q)])
    return {
        "seed":seed,"n_draws":DRAW_COUNT,
        "accepted_burn":burn,
        "burn_attempts":BURN_FACTOR*E,
        "moved_fraction_original_added_edges_after_burn":moved,
        "accepted_swaps_per_spacing_mean":float(np.mean(accepted)),
        "spacing_attempts":SPACING_FACTOR*E,
        "observed_model_analogue_counts":{
            str(q):ess_approx(samples[str(q)]) for q in Q
        },
        "full_draw_history":samples,
    }


def run(source_root):
    model,regions,climate,raw,provenance=source_model(source_root)
    global _LEVEL1_BY_REGION
    _LEVEL1_BY_REGION={i:model["level1"][i] for i in range(365)}
    analog,included=analog_masks(model,regions,climate,raw)
    observed=count_alignment([r["added_bits"] for r in raw],analog)
    assert observed==EXPECTED, f"Original climate source/code mismatch {observed}"
    assert sum(raw[i]["added"] for i in included)==14481
    chains=[run_chain(model,raw,analog,seed) for seed in SEEDS]
    qualified=True
    checks={}
    for q in Q:
        key=str(q)
        a,b=(c["observed_model_analogue_counts"][key] for c in chains)
        agree=abs(a["median"]-b["median"])<=2*max(a["std"],b["std"])
        ess_ok=all(c["observed_model_analogue_counts"][key]["ess_approx"]>=30
                  for c in chains)
        moved_ok=all(c["moved_fraction_original_added_edges_after_burn"]>=.10
                     for c in chains)
        steps_ok=all(c["accepted_swaps_per_spacing_mean"]>=50 for c in chains)
        ok=agree and ess_ok and moved_ok and steps_ok
        qualified &= ok
        checks[key]={"between_chain_median_agreement":agree,
                     "ess_each_ge_30":ess_ok,
                     "moved_fraction_each_ge_0_10":moved_ok,
                     "mean_spacing_accepted_each_ge_50":steps_ok,
                     "qualified":ok}
    comparison={}
    for q in Q:
        key=str(q)
        stats=[c["observed_model_analogue_counts"][key] for c in chains]
        null_med=[x["median"] for x in stats]
        tail=[(sum(v>=observed[key] for v in c["full_draw_history"][key])+1)/100
              for c in chains]
        comparison[key]={
            "observed":observed[key],
            "evaluated_added_butterfly_region_units":14481,
            "null_median_by_chain":null_med,
            "observed_minus_null_median_units_by_chain":[observed[key]-v for v in null_med],
            "observed_minus_null_median_share_pp_by_chain":[100*(observed[key]-v)/14481 for v in null_med],
            "draw_mean_by_chain":[x["mean"] for x in stats],
            "lag1_by_chain":[x["lag1"] for x in stats],
            "ess_by_chain":[x["ess_approx"] for x in stats],
            "nominal_serial_mcmc_right_tail_not_exact_p_by_chain":tail if checks[key]["qualified"] else None,
            "numerically_qualified_exploratory":checks[key]["qualified"]
        }
    return {
        "schema":"chocho_climate_analog_species_region_fixedmargin_null_v0.1",
        "status":"QUALIFIED_EXPLORATORY_MCMC_STRUCTURAL_NULL" if qualified else "NOT_EVALUABLE_SPATIAL_MCMC_MIXING_OR_DIAGNOSTIC",
        "source":"Same pinned author-only WCVP/HOSTS/CHELSA source ZIP, no redistribution",
        "source_sha":provenance["source_frozen_manifest_sha256"],
        "species_in_graph":239,
        "climate_evaluable_species":len(included),
        "source_original_added_units":14553,
        "evaluated_added_units":14481,
        "null_design":"Preserve butterfly×continent and region added-cell degrees, native-host exclusions. Does not preserve introduced plant identities, habitat or butterfly occupancy.",
        "burn_factor":BURN_FACTOR,
        "draws_per_independent_chain":DRAW_COUNT,
        "spacing_factor":SPACING_FACTOR,
        "chains":chains,
        "quality":checks,
        "result":comparison,
        "interpretation_boundary":[
            "Post-hoc after original within-Level1 uniform-reference climate-analogue excess was observed.",
            "This MCMC test conditions on butterfly×continent and region-degree margins, NOT original plant-host geographic networks.",
            "Serial MCMC tails are NOT independent random-permutation p-values or causal introduction effects.",
            "All outcomes are potential butterfly larval host-resource locations, not butterfly occurrence, migration, recruitment or local cofeeding.",
            "No geographic cells' area normalization, historical dates, actual host biomass or observation efforts measured.",
            "The source archive is author-only and its redistribution permissions remain under review."
        ]
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    ans=run(args.source_root)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(ans,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":ans["status"],
        "quality":ans["quality"],
        "result":ans["result"],
        "accepted_burn":[c["accepted_burn"] for c in ans["chains"]],
        "moved_after_burn":[c["moved_fraction_original_added_edges_after_burn"] for c in ans["chains"]],
        "accepted_per_spacing":[c["accepted_swaps_per_spacing_mean"] for c in ans["chains"]]
    },sort_keys=True))

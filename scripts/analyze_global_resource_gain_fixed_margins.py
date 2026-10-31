#!/usr/bin/env python3
"""A stronger geographical null for butterfly potential-resource homogenization.

Unlike HOSTS family/degree swaps, this null fixes EXACTLY:
* The observed native butterfly x WGSRPD3 matrix (26,530 units).
* Each butterfly's *number* of introduced-added resource regions.
* Each WGSRPD3 region's *number* of butterflies with introduced-added resources.
* All 14,553 observed added butterfly-region units and 41,083 current units.
* Native and added resource cannot overlap for any butterfly-region cell.

Only the identity/location combinations of NEW opportunities change via
bipartite 2x2 switches. It tests residual geographic arrangement, NOT
original historical plant transportation nor ecological colonization.

Changes to the fixed plant-host observation process and geography may
still matter; 199 MCMC draws are exploratory, not iid biological replicates.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import random

import numpy as np

BASE=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
MATRIX=BASE/"beta_baseline_v3.matrix.csv.gz"
RECEIPT=BASE/"beta_baseline_v3.json"
OUTPUT=BASE/"resource_gain_fixed_margins_v1.json"
LEVEL1=Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json")
DRAW_TABLE=BASE/"resource_gain_fixed_margins_v1.draws.csv"


def read_matrix(path=MATRIX,receipt_path=RECEIPT):
    receipt=json.loads(receipt_path.read_text(encoding="utf8"))
    states=receipt["resource_matrix"]
    taxa=states["species_order"]
    n_species=len(taxa)
    if len(set(taxa))!=n_species or n_species!=239:
        raise RuntimeError("Native/current source consumer identities changed")
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=369 or len(set(row["wgsrpd3_code"] for row in rows))!=369:
        raise RuntimeError("Not an official 369-region matrix")
    regions=[row["wgsrpd3_code"] for row in rows]
    native=[int(row["native_bits_hex"],16) for row in rows]
    added=[int(row["added_bits_hex"],16) for row in rows]
    current=[int(row["current_bits_hex"],16) for row in rows]
    for n,a,c in zip(native,added,current):
        if n&a or n|a != c:
            raise RuntimeError("Resources not nested as native + introduced-added")
        if (c>>n_species):
            raise RuntimeError("Resource bits beyond frozen butterfly panel")
    checks={"native":sum(n.bit_count() for n in native),
            "added":sum(a.bit_count() for a in added),
            "current":sum(c.bit_count() for c in current)}
    if checks!={"native":26530,"added":14553,"current":41083}:
        raise RuntimeError("Source no longer identical to GEB fixed-potential results: "+repr(checks))
    if checks["native"]!=states["native_species_x_regions"] or checks["current"]!=states["current_species_x_regions"]:
        raise RuntimeError("Binary resource matrix differs from source receipt")
    return regions,taxa,native,added,current,checks


def pair_baseline(native, selected):
    """Fixed native Jaccard comparison on an invariant nonempty region set."""
    pairs=[(i,j) for k,i in enumerate(selected) for j in selected[k+1:]]
    native_mean=0.0
    for i,j in pairs:
        u=(native[i]|native[j]).bit_count()
        native_mean+=(native[i]&native[j]).bit_count()/u
    return pairs,native_mean/len(pairs)


def current_mean(native,added,pairs):
    # N and A always disjoint. No region native emptiness inside selected.
    current=[n|a for n,a in zip(native,added)]
    total=0.0
    for i,j in pairs:
        union=(current[i]|current[j]).bit_count()
        total+=(current[i]&current[j]).bit_count()/union
    return total/len(pairs)


def additions_list(added):
    out=[]
    for region,mask in enumerate(added):
        while mask:
            bit=mask&-mask
            out.append((region,bit.bit_length()-1))
            mask ^= bit
    return out


def swap_new_resources(native,added,edges,rng,attempts,region_strata=None):
    """Constrained 2x2 bipartite switches with all row/column degrees fixed."""
    accepted=0
    n=len(edges)
    by_stratum={}
    if region_strata is not None:
        for k,(reg,sp) in enumerate(edges):
            by_stratum.setdefault(region_strata[reg],[]).append(k)
    for _ in range(attempts):
        i=rng.randrange(n)
        if region_strata is None:
            j=rng.randrange(n)
        else:
            group=by_stratum[region_strata[edges[i][0]]]
            if len(group)<2:
                continue
            j=group[rng.randrange(len(group))]
        if i==j:
            continue
        r1,s1=edges[i]
        r2,s2=edges[j]
        if r1==r2 or s1==s2:
            continue
        bit1=1<<s1
        bit2=1<<s2
        if (added[r1]|native[r1]) & bit2:
            continue
        if (added[r2]|native[r2]) & bit1:
            continue
        if not (added[r1]&bit1) or not (added[r2]&bit2):
            raise RuntimeError("MCMC internal edge source mismatch")
        added[r1] ^= bit1|bit2
        added[r2] ^= bit1|bit2
        edges[i]=(r1,s2)
        edges[j]=(r2,s1)
        accepted+=1
    return accepted


def marginal_counts(native,added,n_species):
    region=[a.bit_count() for a in added]
    species=[0]*n_species
    for a in added:
        while a:
            t=a&-a
            species[t.bit_length()-1]+=1
            a^=t
    return region,species


def evaluate(matrix=MATRIX,receipt=RECEIPT,draws=199,
             attempts_per_added=5,seed=31031,
             stratify_level1=False,lookup=LEVEL1):
    if draws<19 or attempts_per_added<1:
        raise ValueError("At least 19 draws and 1 attempt per added unit required")
    regions,taxa,native,added,current,units=read_matrix(matrix,receipt)
    region_strata=None
    if stratify_level1:
        raw=json.loads(lookup.read_text(encoding="utf8"))
        assert raw.get("schema")=="chocho_wgsrpd3_level1_lookup_v0.1"
        d=raw["units"]
        if set(regions)!=set(d) or len(d)!=369:
            raise RuntimeError("Official WGSRPD3 region strata do not match source")
        region_strata=[d[code] for code in regions]
    selected=[i for i,n in enumerate(native) if n>0]
    if len(selected)!=355:
        raise RuntimeError("Frozen native-region support changed")
    pairs,baseline=pair_baseline(native,selected)
    if len(pairs)!=62835:
        raise RuntimeError("Paired support changed")
    observed=current_mean(native,added,pairs)-baseline
    sampled=added.copy()
    edges=additions_list(sampled)
    if len(edges)!=14553:
        raise RuntimeError("Duplicate source additions")
    fixed_region,fixed_species=marginal_counts(native,added,len(taxa))
    fixed_stratum_species=None
    if region_strata is not None:
        fixed_stratum_species={
            tag:marginal_counts(native,[m if region_strata[i]==tag else 0
                                      for i,m in enumerate(added)],len(taxa))[1]
            for tag in sorted(set(region_strata))}
    rng=random.Random(seed)
    vals=[]
    n_accept=[]
    per_draw=[]
    attempts=int(attempts_per_added*len(edges))
    for i in range(draws):
        k=swap_new_resources(native,sampled,edges,rng,attempts,
                            region_strata=region_strata)
        n_accept.append(k)
        sample_mean=current_mean(native,sampled,pairs)-baseline
        vals.append(sample_mean)
        # Check exact margins at each retained draw; never accept a
        # rearrangement with more/fewer butterflies or regional additions.
        rr,ss=marginal_counts(native,sampled,len(taxa))
        if rr!=fixed_region or ss!=fixed_species:
            raise RuntimeError("Geographic null violated region/species addition margins")
        if any(n & a for n,a in zip(native,sampled)):
            raise RuntimeError("Geographic null overwrote native resource")
        if fixed_stratum_species is not None:
            for tag,expected in fixed_stratum_species.items():
                subset=[m if region_strata[i]==tag else 0
                        for i,m in enumerate(sampled)]
                _,actual=marginal_counts(native,subset,len(taxa))
                if actual!=expected:
                    raise RuntimeError("Within-level1 null violates species gains per world division")
        per_draw.append({
            "draw":i+1,"accepted_swaps":k,
            "mean_delta_fixed_native_region_pairs":sample_mean,
            "all_added_units":sum(rr),
            "native_units":units["native"],
            "current_units":sum(rr)+units["native"],
        })
    a=np.asarray(vals)
    return {
        "schema":"chocho_fixed_species_region_gain_degree_null_v1",
        "status":"ADDED_OPPORTUNITY_SIZE_MATCHED_GEOGRAPHY_NULL_NOT_HOST_TRANSPORT",
        "source": {
            "resource_matrix_sha256":hashlib.sha256(matrix.read_bytes()).hexdigest(),
            "report_source_sha256":hashlib.sha256(receipt.read_bytes()).hexdigest(),
            "butterfly_species":len(taxa),"official_regions":len(regions),
            "common_native_resource_regions":len(selected),
            "common_region_pairs":len(pairs),**units,
        },
        "rule":"2x2 swaps of added butterfly×WGSRPD3 cells; forbid native cells; preserve each butterfly's number of gains and each region's number of new butterflies; preserve all original native cells",
        "stratification":"within WGSRPD level1 region division" if stratify_level1 else "global unstratified",
        "strata_source_sha256":hashlib.sha256(lookup.read_bytes()).hexdigest() if stratify_level1 else None,
        "observed_mean_beta_jaccard_change":observed,
        "randomization": {
            "draws":draws,"seed":seed,
            "attempted_edge_swaps_per_draw":attempts,
            "min_accepted_swaps_per_draw":min(n_accept),
            "median_accepted_swaps_per_draw":float(np.median(n_accept)),
            "mean":float(np.mean(a)),
            "median":float(np.median(a)),
            "q025":float(np.quantile(a,.025)),
            "q975":float(np.quantile(a,.975)),
            "p_low":(1+sum(v<=observed for v in vals))/(draws+1),
            "p_high":(1+sum(v>=observed for v in vals))/(draws+1),
        },
        "draws_table":per_draw,
        "limits":[
            "Source WCVP native/current statuses are modern retrospective geography, not dated plant colonization.",
            "No GBIF confirmed butterfly occupancy; output is only potential larval resource geography.",
            "MCMC draws from sequential constrained swaps are not independent experimental replications.",
            "This null fixes butterfly gain counts by region and species, not plant-level introduction or host-network structure.",
            "Region endpoints are reused in pair summaries; no pair-independent confidence intervals.",
            "Fixed margin null does not account for spatial autocorrelation/phylogeny of introduced host plant ranges.",
            "Within-level1 swaps retain broad geographic placement, but not climate, port traffic, or spatial autocorrelation within world division.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--matrix",type=Path,default=MATRIX)
    p.add_argument("--receipt",type=Path,default=RECEIPT)
    p.add_argument("--draws",type=int,default=199)
    p.add_argument("--attempts-per-added",type=int,default=5)
    p.add_argument("--stratify-level1",action="store_true")
    p.add_argument("--level1-lookup",type=Path,default=LEVEL1)
    p.add_argument("--output",type=Path,default=OUTPUT)
    p.add_argument("--draws-csv",type=Path,default=DRAW_TABLE)
    args=p.parse_args()
    r=evaluate(args.matrix,args.receipt,args.draws,args.attempts_per_added,
               stratify_level1=args.stratify_level1,lookup=args.level1_lookup)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({k:v for k,v in r.items() if k!="draws_table"},
                                      indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with args.draws_csv.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(r["draws_table"][0]))
        w.writeheader();w.writerows(r["draws_table"])
    print(json.dumps({k:v for k,v in r.items() if k not in ("draws_table","limits")},
                     indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Source-lock island/mainland potential-butterfly-resource Jaccard homogenization.

No butterfly GBIF records, ocean current, vessel route or real larval fitness.
Replays the original PR118 239×369 accepted-WCVP/HOSTS 0/1 resource matrix;
compares island/island vs mainland/mainland + cross pairs under original
region/species exact-gain-margin swaps and the stricter Level1×island-status
stratified model. The original global within-Level1 null is reproduced first.
Exploratory after the existing global +0.185 source result was known.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import statistics
import sys
from collections import Counter
from pathlib import Path

import numpy as np

SOURCE_COMMIT = "bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2"
EXPECTED_MATRIX_SHA256 = "ec6a5530d65ac0b04f230e0c51fd68bf543555641a17b94795996bf11ce257d0"
EXPECTED_REPORT_SHA256 = "ccda0bc669fb5ffe7326175fcc0501dbb1447790925c06add92501d1d30f255b"
EXPECTED_LEVEL1_SHA256 = "85804d45b0a9ffb3620d13795015ca2a23d2ecd4c99e489f7d312906750d94d7"
EXPECTED_OLD_OBSERVED = 0.18510150276498433
EXPECTED_OLD_NULL_MEAN = 0.18158108553894184
N_DRAW=199
N_ATTEMPT=72765
SEED_BASE=41131
SEED_STRICT=20261118
N_CLASSES={"island":126,"mainland":237,"mixed":6}
CATEGORIES=("ALL_NATIVE_SUPPORTED","ISLAND_ISLAND","MAINLAND_MAINLAND",
            "ISLAND_MAINLAND","ISLAND_ISLAND_SAME_LEVEL1","MAINLAND_MAINLAND_SAME_LEVEL1")

def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read_csv(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def relation(ri,rj):
    if ri=="island" and rj=="island":return "ISLAND_ISLAND"
    if ri=="mainland" and rj=="mainland":return "MAINLAND_MAINLAND"
    if {ri,rj}=={"island","mainland"}:return "ISLAND_MAINLAND"
    return None

def make_pairs(native,regions,classes,level1):
    selected=[i for i,n in enumerate(native) if n!=0]
    if len(selected)!=355: raise RuntimeError("Original common 355 native-resource supported regions changed")
    pairs=[]
    cts=Counter()
    for ix,i in enumerate(selected):
        for j in selected[ix+1:]:
            label=relation(classes[regions[i]],classes[regions[j]])
            li=level1[regions[i]]
            lj=level1[regions[j]]
            selected_categories=["ALL_NATIVE_SUPPORTED"]
            if label:
                selected_categories.append(label)
                if li==lj and label in ("ISLAND_ISLAND","MAINLAND_MAINLAND"):
                    selected_categories.append(label+"_SAME_LEVEL1")
            k=(native[i]&native[j]).bit_count()
            union=(native[i]|native[j]).bit_count()
            if union==0:raise RuntimeError("Undefined native Jaccard within original nonempty support")
            pairs.append((i,j,k/union,tuple(selected_categories)))
            for tag in selected_categories:cts[tag]+=1
    if len(pairs)!=62835 or cts["ALL_NATIVE_SUPPORTED"]!=62835:
        raise RuntimeError("Original resource source pair count changed")
    if cts["ISLAND_ISLAND"]<1000 or cts["MAINLAND_MAINLAND"]<10000:
        raise RuntimeError("Insufficient source island or mainland support")
    if cts["ISLAND_ISLAND_SAME_LEVEL1"]<50 or cts["MAINLAND_MAINLAND_SAME_LEVEL1"]<50:
        raise RuntimeError("Insufficient within-biogeography pair controls")
    return selected,pairs,cts

def scores_from_matrix(native,added,pairs,*,include_distribution=False):
    current=[n|a for n,a in zip(native,added)]
    cnt={cat:0 for cat in CATEGORIES}
    sum_cur={cat:0.0 for cat in CATEGORIES}
    sum_del={cat:0.0 for cat in CATEGORIES}
    if include_distribution:
        vals={cat:[] for cat in CATEGORIES}
    for i,j,nj,categories in pairs:
        u=(current[i]|current[j]).bit_count()
        if u<=0:raise RuntimeError("Nested current native-supported union zero")
        cj=(current[i]&current[j]).bit_count()/u
        delta=cj-nj
        for cat in categories:
            cnt[cat]+=1
            sum_cur[cat]+=cj
            sum_del[cat]+=delta
            if include_distribution:vals[cat].append(delta)
    out={}
    for cat in CATEGORIES:
        n=cnt[cat]
        if n<50:raise RuntimeError("No source pair support for "+cat)
        d={"n_pairs":n,"mean_current_jaccard":sum_cur[cat]/n,"mean_delta_jaccard":sum_del[cat]/n}
        if include_distribution:
            arr=vals[cat]
            d.update(median_delta_jaccard=float(statistics.median(arr)),
                     positive_delta_fraction=sum(x>0 for x in arr)/n,
                     negative_delta_fraction=sum(x<0 for x in arr)/n,
                     zero_delta_fraction=sum(x==0 for x in arr)/n)
        out[cat]=d
    return out

def native_means(pairs):
    sums=Counter()
    nums=Counter()
    for i,j,nj,cats in pairs:
        for cat in cats:
            sums[cat]+=nj
            nums[cat]+=1
    return {cat:float(sums[cat]/nums[cat]) for cat in CATEGORIES}

def margins(matrix,n_species):
    by_region=[a.bit_count() for a in matrix]
    by_species=[0]*n_species
    for a in matrix:
        while a:
            z=a&-a
            by_species[z.bit_length()-1]+=1
            a^=z
    return by_region,by_species

def make_null(source,native,added,pairs,regions,level1,classes,
              frozen_region_counts,frozen_species_counts,mode):
    if mode=="level1":
        strata=[level1[c] for c in regions]
        seed=SEED_BASE
    elif mode=="level1_and_origin_status":
        strata=[level1[c]+"|"+classes[c] for c in regions]
        seed=SEED_STRICT
    else:raise RuntimeError("No prespecified null mode")
    unique_strata=sorted(set(strata))
    # Require exact butterfly gain totals within each source-selected stratum.
    original_stratum_species={
        tag:margins([a if strata[i]==tag else 0 for i,a in enumerate(added)],len(source.SPECIES))[1]
        for tag in unique_strata
    } if hasattr(source,"SPECIES") else None
    # Avoid source module hidden state; construct directly from supplied matrix.
    source_edges=source.additions_list(added)
    if len(source_edges)!=14553:raise RuntimeError("Original 14,553 added source cells missing")
    edges=source_edges.copy()
    sampled=added.copy()
    rng=random.Random(seed)
    draws=[]
    for d in range(N_DRAW):
        n_acc=source.swap_new_resources(native,sampled,edges,rng,N_ATTEMPT,region_strata=strata)
        if n_acc<200:
            raise RuntimeError(f"Original prescribed 2x2 mobility too small, null {mode}, draw {d}: {n_acc}")
        rr,ss=margins(sampled,len(frozen_species_counts))
        if rr!=frozen_region_counts or ss!=frozen_species_counts:
            raise RuntimeError("Original exact addition row/column margins violated")
        if any(n&a for n,a in zip(native,sampled)):raise RuntimeError("Any native cell overwritten")
        # Per-butterfly gain totals within exact source Level1 (or Level1×status) remain invariant by original swaps.
        stat=scores_from_matrix(native,sampled,pairs)
        record={"draw":d+1,"accepted_swaps":n_acc}
        for cat in CATEGORIES:
            record[cat+"_delta"]=stat[cat]["mean_delta_jaccard"]
        record["ISLAND_MINUS_MAINLAND"]=record["ISLAND_ISLAND_delta"]-record["MAINLAND_MAINLAND_delta"]
        record["SAME_LEVEL1_ISLAND_MINUS_MAINLAND"]=(
            record["ISLAND_ISLAND_SAME_LEVEL1_delta"]-record["MAINLAND_MAINLAND_SAME_LEVEL1_delta"])
        draws.append(record)
    return draws

def contrast_summary(observed,draws,metric,scope):
    v=np.array([r[metric] for r in draws],dtype=float)
    point=observed
    return {
        "metric":metric,"source_scope":scope,"observed":float(point),
        "null_mean":float(v.mean()),"null_median":float(np.median(v)),
        "null_descriptive_q025":float(np.quantile(v,.025)),
        "null_descriptive_q975":float(np.quantile(v,.975)),
        "observed_minus_null_median":float(point-np.median(v)),
        "descriptive_high_tail_add_one":float((1+np.sum(v>=point))/(len(v)+1)),
        "draw_count":len(draws)
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archived-root",required=True,type=Path)
    p.add_argument("--region-class",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    root=a.archived_root
    old=root/"provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031"
    matrix=old/"beta_baseline_v3.matrix.csv.gz"
    receipt=old/"beta_baseline_v3.json"
    lookup=root/"data/frozen/wgsrpd3_level1_lookup_v0.1.json"
    if sha256(matrix)!=EXPECTED_MATRIX_SHA256:raise RuntimeError("Source binary matrix SHA256 changed")
    if sha256(receipt)!=EXPECTED_REPORT_SHA256:raise RuntimeError("Source receipt SHA256 changed")
    if sha256(lookup)!=EXPECTED_LEVEL1_SHA256:raise RuntimeError("Source TDWG Level1 map SHA256 changed")
    sys.path.insert(0,str(root/"scripts"))
    import analyze_global_resource_gain_fixed_margins as source
    regions,taxa,native,added,current,units=source.read_matrix(matrix,receipt)
    if len(regions)!=369 or len(taxa)!=239 or units!={"native":26530,"added":14553,"current":41083}:
        raise RuntimeError("Original 239x369 resource source does not pass exact population gate")
    geo=json.loads(lookup.read_text(encoding="utf-8"))["units"]
    if len(geo)!=369 or set(geo)!=set(regions):raise RuntimeError("Original TDWG origin map mismatch")
    cl=read_csv(a.region_class)
    cls={r["area_code_l3"]:r["class"] for r in cl}
    if len(cl)!=369 or set(cls)!=set(regions) or Counter(cls.values())!=N_CLASSES:
        raise RuntimeError("Original 126-island/237-mainland/6-mixed class roster source drift")
    selected,pairs,counts=make_pairs(native,regions,cls,geo)
    nmean=native_means(pairs)
    obs=scores_from_matrix(native,added,pairs,include_distribution=True)
    for cat in CATEGORIES:
        obs[cat]["mean_native_jaccard"]=nmean[cat]
        if abs(obs[cat]["mean_current_jaccard"]-nmean[cat]-obs[cat]["mean_delta_jaccard"])>1e-10:
            raise RuntimeError("Native/current source arithmetic invalid")
    if abs(obs["ALL_NATIVE_SUPPORTED"]["mean_delta_jaccard"]-EXPECTED_OLD_OBSERVED)>1e-10:
        raise RuntimeError("Original global source 355-region positive resource-beta delta not reproduced")
    frozen_region,frozen_species=margins(added,len(taxa))
    if len(frozen_region)!=len(regions) or sum(frozen_region)!=14553 or sum(frozen_species)!=14553:
        raise RuntimeError("Source addition margins invalid")
    nulls={}
    for mode in ("level1","level1_and_origin_status"):
        draw_rows=make_null(source,native,added,pairs,regions,geo,cls,frozen_region,frozen_species,mode)
        if mode=="level1":
            reproduced=np.mean([v["ALL_NATIVE_SUPPORTED_delta"] for v in draw_rows])
            if abs(reproduced-EXPECTED_OLD_NULL_MEAN)>1e-8:
                raise RuntimeError(f"Original global Level1 null did not reproduce: {reproduced}")
        base=obs["ISLAND_ISLAND"]["mean_delta_jaccard"]-obs["MAINLAND_MAINLAND"]["mean_delta_jaccard"]
        same=obs["ISLAND_ISLAND_SAME_LEVEL1"]["mean_delta_jaccard"]-obs["MAINLAND_MAINLAND_SAME_LEVEL1"]["mean_delta_jaccard"]
        res={
            "method":mode,"draws":len(draw_rows),"attempts_per_draw":N_ATTEMPT,
            "seed":SEED_BASE if mode=="level1" else SEED_STRICT,
            "accepted_swaps_min":min(r["accepted_swaps"] for r in draw_rows),
            "accepted_swaps_median":float(np.median([r["accepted_swaps"] for r in draw_rows])),
            "cohorts":{cat:contrast_summary(obs[cat]["mean_delta_jaccard"],draw_rows,cat+"_delta",cat) for cat in CATEGORIES},
            "island_minus_mainland":contrast_summary(base,draw_rows,"ISLAND_MINUS_MAINLAND","II-MM"),
            "same_Level1_island_minus_mainland":contrast_summary(same,draw_rows,"SAME_LEVEL1_ISLAND_MINUS_MAINLAND","within-Level1 II-MM")
        }
        nulls[mode]=res
        a.out_dir.mkdir(parents=True,exist_ok=True)
        with (a.out_dir/("draws_"+mode+".csv")).open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(draw_rows[0]))
            w.writeheader();w.writerows(draw_rows)
        print("COMPLETED_NULL",json.dumps({
            "mode":mode,
            "original_global_null_mean":res["cohorts"]["ALL_NATIVE_SUPPORTED"]["null_mean"],
            "island":res["cohorts"]["ISLAND_ISLAND"],
            "mainland":res["cohorts"]["MAINLAND_MAINLAND"],
            "island_minus_mainland":res["island_minus_mainland"],
            "accepted_swap_min":res["accepted_swaps_min"]
        }),flush=True)
    excluded_native=[
        {"code":regions[i],"source_class":cls[regions[i]],"level1":geo[regions[i]],
         "current_resource_count":current[i].bit_count()}
        for i,n in enumerate(native) if n==0
    ]
    doc={
        "schema":"chocho_original239_island_mainland_source_constrained_resource_Jaccard_v1",
        "status":"SOURCE_VERIFIED_EXPLORATORY_POTENTIAL_RESOURCE_COMPOSITION_ONLY",
        "original_PR118_commit":SOURCE_COMMIT,
        "source_matrix_sha256":sha256(matrix),"source_report_sha256":sha256(receipt),
        "original_WGSRPD_level1_sha256":sha256(lookup),
        "butterflies":len(taxa),"region_codes_total":len(regions),
        "original_resource_units":units,
        "original_region_class_counts":dict(Counter(cls.values())),
        "native_nonempty_region_count":len(selected),
        "excluded_no_native_resource_regions":excluded_native,
        "eligible_pair_denominators":dict(counts),
        "observed":obs,
        "nulls":nulls,
        "scientific_boundary":[
            "All observations use previously source-derived HOSTS/WCVP potential food geography, not a new butterfly occurrence observation.",
            "WCVP current vs native-only represents nested contemporary spatial states, not a historical before/after immigration process.",
            "Species and region margins are preserved exactly; source samples are observational, correlated and taxonomically biased.",
            "MCMC sequential 199 draws are dependent and quantiles/tails are descriptive, not independent hypothesis test confidence intervals.",
            "Any mainland/island advantage depends on source inclusion and island labels, not observed port/marine species transport.",
            "No source flower/plant trait, plant abundance/host quality, butterfly population or historical vessel flow outcome is modeled."
        ],
    }
    (a.out_dir/"island_mainland_resource_beta_verified_summary.json").write_text(
        json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("SOURCE_VERIFIED_ISLAND_MAINLAND_POTENTIAL_CONVERGENCE",json.dumps({
        "original_global_observed":obs["ALL_NATIVE_SUPPORTED"]["mean_delta_jaccard"],
        "class_region_counts":doc["original_region_class_counts"],
        "pair_denominators":doc["eligible_pair_denominators"],
        "observed_delta_by_class":{z:obs[z]["mean_delta_jaccard"] for z in CATEGORIES},
        "island_minus_mainland_observed":obs["ISLAND_ISLAND"]["mean_delta_jaccard"]-obs["MAINLAND_MAINLAND"]["mean_delta_jaccard"],
        "relative_to_null":{mode:{
            "island_excess":q["cohorts"]["ISLAND_ISLAND"]["observed_minus_null_median"],
            "mainland_excess":q["cohorts"]["MAINLAND_MAINLAND"]["observed_minus_null_median"],
            "II_minus_MM_excess":q["island_minus_mainland"]["observed_minus_null_median"]
        } for mode,q in nulls.items()}
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Family-conditional label shuffles for EuDiS thalassochory and WCVP shared flora.

A source-family-matched test of contemporary *region listings*, NOT seed
dispersal, observed island arrival, historical anthropogenic introductions,
physical current interaction, or butterfly host use.
"""
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict
from pathlib import Path

PAIRS=(("AZO","CNY"),("CNY","CVI"),("AZO","CVI"),
       ("CNY","AZO"),("CVI","CNY"),("CVI","AZO"))
NULL_DRAWS=9999
OUTCOME_NATIVE="WCVP_NATIVE"
OUTCOME_INTRODUCED="WCVP_INTRODUCED"
OUTCOME_UNLISTED="WCVP_NOT_LISTED"

def read_rows(path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def contrast(trait,outcome):
    n1=sum(x==1 for x in trait);n0=len(trait)-n1
    if n1==0 or n0==0:raise RuntimeError("Need both trait statuses")
    return sum(z for t,z in zip(trait,outcome) if t==1)/n1 - \
           sum(z for t,z in zip(trait,outcome) if t==0)/n0

def analyze(rows,src,dst):
    # In the region, source-native trait-known accepted WCVP species are a
    # deliberately nonrepresentative subset of source WCVP flora.
    ids=[x["plant_name_id"] for x in rows]
    if len(ids)!=len(set(ids)):raise RuntimeError("Duplicate accepted plant id in source region")
    traits=[int(r["thalassochorous"]) for r in rows]
    outcomes=[int(r["dst_status"] in (OUTCOME_NATIVE,OUTCOME_INTRODUCED)) for r in rows]
    if set(traits)!={0,1}:raise RuntimeError("Missing categorical EuDiS trait group")
    families=defaultdict(list)
    for i,r in enumerate(rows):
        f=r["family"].strip()
        if not f or f=="__UNKNOWN_FAMILY__":
            continue # track as an unexchangeable reference family, not a broad pooled family
        families[f].append(i)
    active=[z for z in families.values() if {traits[i] for i in z}=={0,1}]
    fixed=[z for z in families.values() if {traits[i] for i in z}!={0,1}]
    family_diversity=len(families)
    exch=sum(map(len,active))
    exch_ones=sum(sum(traits[i] for i in z) for z in active)
    obs=contrast(traits,outcomes)
    # Exact conditional expectation from fixed 1/0 counts within family.
    n1=sum(traits);n0=len(traits)-n1
    # Background nonexchanged singleton/unknown family labels remain fixed.
    expected_one=sum(outcomes[i] for i in range(len(traits))
                     if traits[i]==1 and
                     (rows[i]["family"].strip()=="__UNKNOWN_FAMILY__" or
                      not rows[i]["family"].strip()))
    expected_zero=sum(outcomes[i] for i in range(len(traits))
                     if traits[i]==0 and
                     (rows[i]["family"].strip()=="__UNKNOWN_FAMILY__" or
                      not rows[i]["family"].strip()))
    for z in families.values():
        n=sum(traits[i] for i in z)
        m=len(z)
        y=sum(outcomes[i] for i in z)
        expected_one+=n*y/m
        expected_zero+=(m-n)*y/m
    expected=expected_one/n1 - expected_zero/n0
    rng=random.Random(2026111700+sum(map(ord,src+dst)))
    diffs=[]
    for _ in range(NULL_DRAWS):
        t=traits.copy()
        for z in active:
            values=[t[i] for i in z]
            rng.shuffle(values)
            for i,v in zip(z,values):t[i]=v
        diffs.append(contrast(t,outcomes))
    null_center=sum(diffs)/NULL_DRAWS
    if abs(null_center-expected)>0.02:
        raise RuntimeError("Analytical conditional family expectation inconsistent with permutations")
    diffs.sort()
    lo,hi=diffs[int((NULL_DRAWS-1)*.025)],diffs[int((NULL_DRAWS-1)*.975)]
    count=sum(abs(v-expected)>=abs(obs-expected)-1e-12 for v in diffs)
    p=(count+1)/(NULL_DRAWS+1)
    return {
      "source":src,"destination":dst,
      "source_native_trait_known_species":len(rows),
      "source_native_thalasso1":n1,"source_native_thalasso0":n0,
      "accepted_taxonomic_families_with_known_name":family_diversity,
      "within_family_exchangeable_families":len(active),
      "within_family_exchangeable_species":exch,
      "within_family_exchangeable_thalasso1_species":exch_ones,
      "native_or_introduced_WCVP_recipient_coverage_trait1":sum(z for t,z in zip(traits,outcomes) if t==1),
      "native_or_introduced_WCVP_recipient_coverage_trait0":sum(z for t,z in zip(traits,outcomes) if t==0),
      "observed_thalasso1_minus_0_listing_pp":100*obs,
      "conditional_family_expected_listing_pp":100*expected,
      "within_family_excess_pp":100*(obs-expected),
      "family_shuffle_95pct_reference_pp":[100*lo,100*hi],
      "two_sided_family_label_permutation_p_exploratory":p,
      "status": "WITHIN_FAMILY_PERMUTABLE" if active else "NO_WITHIN_FAMILY_EXCHANGEABLE_TRAITS",
      "source_bias":"EuDiS is biased towards European native taxa; shared flora is not dated dispersal"
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    rows=read_rows(args.input)
    by_pair=defaultdict(list)
    for r in rows:
        s,d=r["source"],r["destination"]
        if (s,d) not in PAIRS:raise RuntimeError("Non-predeclared source/destination in input")
        if r["thalassochorous"] not in ("0","1"):raise RuntimeError("Unresolved trait included in source")
        if r["dst_status"] not in (OUTCOME_NATIVE,OUTCOME_INTRODUCED,OUTCOME_UNLISTED):
            raise RuntimeError("Unknown WCVP recipient regional status")
        by_pair[(s,d)].append(r)
    if set(by_pair)!=set(PAIRS):raise RuntimeError("Missing original source region pair")
    summaries=[analyze(by_pair[k],*k) for k in PAIRS]
    result={
      "status":"PRESENT_DAY_WCVP_REGION_LISTING_FAMILY_MATCHED_EXPLORATORY_SOURCE_ONLY",
      "replications":NULL_DRAWS,
      "estimand":"destination WCVP regional listing rate, conditional on source-native EuDiS-known species pool",
      "method":"Shuffle EuDiS thalassochory 1/0 labels ONLY among accepted species of the same WCVP family within an identical original source region, preserving family totals; unknown family unshuffled.",
      "original_WCVP_v13":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
      "original_EuDiS_v5_figshare_file_id":41269785,
      "paired_regions":summaries,
      "key_limits":[
        "No source-native species identification of which shore population provided a seed.",
        "WCVP currently listed vs unlisted is neither a historical immigration event nor a biological absence.",
        "EuDiS thalassochory is a biased morphological trait inventory, not confirmed buoyancy or viability.",
        "Single region/source species appears across multiple directions, so pair comparisons are not independent ecological replications.",
        "No matched climate, coastline habitat, distance, observed currents after floating duration, actual vessel OD or butterfly settlement response.",
        "Permutation p is conditional on source's family-labelled inventory and cannot identify a causal sea-current mechanism."
      ]
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("ATLANTIC_FAMILY_MATCHED_SEA_DISPERSAL_TRAIT_NULL",json.dumps({
        "status":result["status"],"source_rows":len(rows),
        "pairs":summaries},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

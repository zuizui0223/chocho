#!/usr/bin/env python3
"""Exploratory plant-family-matched comparison of EuDiS thalassochory
versus WCVP destination-native reporting. Retains exact source 3 pair scopes.
Not a sea-current, seed-tracking, ship-transport or butterfly colonization test.
"""
from __future__ import annotations
import argparse
import csv
import json
import random
import statistics
from collections import defaultdict,Counter
from pathlib import Path

ROUTES=(("AZO","CNY"),("CNY","CVI"),("AZO","CVI"))
SEED=20261117

def load(path):
    with path.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def percentile(samples,p):
    z=sorted(samples)
    k=(len(z)-1)*p
    lo=int(k);hi=min(lo+1,len(z)-1)
    return z[lo]+(z[hi]-z[lo])*(k-lo)

def estimate(families):
    # Family-matched excess, standardized to the thalassochorous source
    # family's distribution, not a pooled between-family ecological signal.
    n1=0;total=0.0
    for f,rs in families.items():
        a=[r["native"] for r in rs if r["trait"]==1]
        b=[r["native"] for r in rs if r["trait"]==0]
        if not (a and b):continue
        n1+=len(a)
        total+=len(a)*(statistics.mean(a)-statistics.mean(b))
    return total/n1 if n1 else None

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-csv",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    rows=load(a.source_csv)
    result={"status":"EXPLORATORY_FAMILY_MATCHED_NATIVE_PLANT_RECORDING_NOT_SEED_TRANSPORT",
      "source_CI":43898725162,"source_species_cohorts":"Only WCVP accepted source-native species with EuDiS v5 identified 0/1",
      "route_scope":"AZO->CNY CNY->CVI AZO->CVI fixed",
      "permutation_seed":SEED,"permutations":4999,
      "end_point":"WCVP destination-native catalogue record (not arrival history)",
      "results":[],"caveats":[
        "Only families with both thalassochory 0 and 1 have within-family identification; unmatched positive families are not silently assigned 0.",
        "Permutation is a sensitivity analysis after looking at source sign; not a pre-registered confirmatory test.",
        "Shared native flora and EuDiS category are observational; climatic and geographic filtering and phylogeny within family remain unadjusted.",
        "The same plant species may recur across route comparisons, so tests are not independent.",
        "Unreported WCVP species at destination are classified as not indexed, NOT proven biological absence.",
        "No physical water-drift, ship OD, introduction chronology, butterfly use, or colonization outcome tested.",
      ]}
    for src,dst in ROUTES:
        sub=[r for r in rows if r["source"]==src and r["destination"]==dst]
        if not sub:raise RuntimeError(f"No pinned source for {src}->{dst}")
        fam=defaultdict(list)
        tr_count=Counter()
        for r in sub:
            tr=int(r["EuDiS_thalassochorous"])
            y=int(r["destination_origin"]=="destination_native")
            if tr not in (0,1):raise RuntimeError("Source trait unrecognized")
            family=str(r["family"]).strip() or "__unknown__"
            fam[family].append({"trait":tr,"native":y})
            tr_count[tr]+=1
        valid={k:v for k,v in fam.items() if any(x["trait"]==0 for x in v) and any(x["trait"]==1 for x in v)}
        raw_a=[r for rs in fam.values() for r in rs if r["trait"]==1]
        raw_b=[r for rs in fam.values() for r in rs if r["trait"]==0]
        if (len(raw_a),len(raw_b)) not in {(24,354),(51,608)}:
            raise RuntimeError("Original source-positive/negative counts changed")
        observed=estimate(valid)
        pos=sum(sum(x["trait"]==1 for x in v) for v in valid.values())
        neg=sum(sum(x["trait"]==0 for x in v) for v in valid.values())
        eligible=len(valid)>=3 and pos>=10 and neg>=10
        rec={
          "source":src,"destination":dst,
          "source_positive_species":len(raw_a),"source_negative_species":len(raw_b),
          "raw_native_fraction_positive":statistics.mean(r["native"] for r in raw_a),
          "raw_native_fraction_negative":statistics.mean(r["native"] for r in raw_b),
          "raw_native_sharing_difference_pp":100*(statistics.mean(r["native"] for r in raw_a)-
                                                  statistics.mean(r["native"] for r in raw_b)),
          "source_family_count":len(fam),"both_trait_state_family_count":len(valid),
          "family_match_positive_species":pos,"family_match_negative_species":neg,
          "family_matched_native_fraction_difference_pp":100*observed if observed is not None else None,
          "statistical_gate":"EXPLORATORY_PERMUTATION_REPORTED" if eligible else "INSUFFICIENT_WITHIN_FAMILY_SUPPORT",
          "not_independent_transport_test":True,
        }
        if eligible:
            rng=random.Random(SEED)
            values=[]
            for _ in range(4999):
                res={}
                for f,grp in valid.items():
                    y=[r["native"] for r in grp]
                    rng.shuffle(y)
                    res[f]=[{"trait":r["trait"],"native":z} for r,z in zip(grp,y)]
                values.append(estimate(res))
            assert observed is not None
            rec["family_permutation_null_95pp"]=[100*percentile(values,.025),100*percentile(values,.975)]
            rec["family_permutation_two_sided_exploratory_p"]=(1+sum(abs(v)>=abs(observed)-1e-12 for v in values))/(len(values)+1)
        result["results"].append(rec)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ATLANTIC_WCVP_FAMILY_MATCHED_THALASSO_TEST",json.dumps(result["results"],ensure_ascii=False),flush=True)

if __name__=="__main__":main()

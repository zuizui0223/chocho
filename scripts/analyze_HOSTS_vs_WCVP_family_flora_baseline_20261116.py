#!/usr/bin/env python3
"""Is excess HOSTS alien-plant reporting specific to host families or just flora?

WCVP's accepted contemporary wild vascular flora is a geographic taxon pool,
NOT actual larval food abundance, a random encounter model or a field experiment.
All results are exploratory and must retain this distinction.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.stats import mannwhitneyu
from analyze_HOSTS_reported_use_island_fao_ports_20261116 import build_rows, unweighted_lm

MODES=("all_hosts","exclude_FAO_exact","exclude_FAO_exact_plus_genus_wildcards",
       "exclude_FAO_all_crop_genera")
COHORTS=("HOSTS_Papilionoidea","fixed239_butterflies","all_HOSTS_lepidoptera")


def read(p):
    with open(p,newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))


def check_flora(path, island_csv):
    # WCVP full accepted flora counts must reproduce pre-existing published
    # source-derived island flora counts; fail closed on origin/taxon drift.
    f=defaultdict(lambda:{"native":0,"introduced":0})
    species=defaultdict(lambda:{"native":0,"introduced":0})
    for z in read(path):
        if z["cropmode"]!="all_hosts":
            continue
        flag="introduced" if z["introduced"]=="1" else "native"
        if z["family"]=="__ALL_FAMILIES__":
            f[z["area_code_l3"]][flag]+=int(z["n"])
        else:
            species[(z["area_code_l3"],z["family"])][flag]+=int(z["n"])
    old={r["code"]:r for r in read(island_csv)}
    both=[c for c in f if c in old]
    if len(both)<35:raise ValueError("Too few original island-flora overlap controls")
    mismatch=[]
    for c in both:
        for flag in ("native","introduced"):
            if f[c][flag]!=int(old[c][flag]):
                mismatch.append([c,flag,f[c][flag],old[c][flag]])
    if mismatch:raise RuntimeError(f"Accepted WCVP source flora changed: {mismatch[:10]}")
    return {"matching_island_regions":len(both),"mismatches":0,
            "all_source_baseline_regions":len(f)}


def baseline(path):
    out=defaultdict(lambda:[0,0])
    for x in read(path):
        k=(x["area_code_l3"],x["cropmode"],x["family"])
        out[k][int(x["introduced"])]+=int(x["n"])
    return out


def select(rows,cohort,mode):
    for r in rows:
        if r["mapping_tier"]!="exact" or "|" in r["codes"]:
            continue  # only exact single-WGSRPD3 geography, never composite
        if r["rearing"]=="Y":
            continue
        if cohort=="HOSTS_Papilionoidea" and not r["papilionoidea"]:
            continue
        if cohort=="fixed239_butterflies" and not r["focal_butterfly239"]:
            continue
        if mode!="all_hosts" and r["crop_exact"]:
            continue
        if mode in MODES[2:] and r["crop_genus_wildcard"]:
            continue
        if mode=="exclude_FAO_all_crop_genera" and r["crop_any_FAO_genus"]:
            continue
        yield r


def sample(rows,flora,name_resolution,minimum,mode):
    fam_by_id=defaultdict(set)
    for r in read(name_resolution):
        if r["match_status"]=="matched":
            fam_by_id[r["accepted_plant_name_id"]].add(r["family"])
    grouped=defaultdict(list)
    for r in rows:grouped[r["location"]].append(r)
    data=[]
    excluded=Counter()
    for loc,arr in sorted(grouped.items()):
        if len(arr)<minimum:continue
        original_n=len(arr)
        codes={r["codes"] for r in arr}
        if len(codes)!=1 or "|" in next(iter(codes)):
            raise ValueError("Exact one-source geography promised but violated: "+loc)
        code=next(iter(codes))
        allv=flora.get((code,mode,"__ALL_FAMILIES__"),[0,0])
        if sum(allv)==0:excluded["no_flora_baseline"]+=1;continue
        use=[]
        for r in arr:
            f=fam_by_id.get(r["accepted_plant_id"],set())
            if len(f)!=1:excluded["unknown_or_ambiguous_host_family"]+=1;continue
            fam=next(iter(f))
            fcounts=flora.get((code,mode,fam),[0,0])
            if sum(fcounts)==0:
                excluded["no_accepted_family_flora"]+=1
                continue
            use.append((r,int(r["origin"]=="introduced"),fcounts[1]/sum(fcounts)))
        if len(use)<minimum:
            excluded["insufficient_source_rows_family_matched"]+=1
            continue
        n=len(use);i=sum(v for _,v,_ in use)
        expected=float(np.mean([p for _,_,p in use]))
        total_share=allv[1]/sum(allv)
        d={
            "location":loc,"WGSRPD3_code":code,
            "origin_class":arr[0]["origin_class"],
            "is_island":int(arr[0]["comparison_group"]=="island"),
            "WPI_ports":int(arr[0]["WPI_ports"]),
            "region_area_km2":float(arr[0]["region_area_km2"]),
            "level1":arr[0]["level1"],
            "n":n,"original_location_rows_before_family_match":original_n,
            "introduced_host_reports":i,"native_host_reports":n-i,
            "host_report_alien_share":i/n,
            "flora_all_alien_share":total_share,
            "flora_family_mixture_expected_share":expected,
            "host_minus_flora_all":i/n-total_share,
            "host_minus_flora_family":i/n-expected,
            "log_odds_host_minus_flora_all":
                math.log((i+.5)/(n-i+.5))-math.log((allv[1]+.5)/(allv[0]+.5)),
            "total_record_count":n,
        }
        data.append(d)
    return data,dict(excluded)


def quantile_ci(vals,rng,N=4999):
    if len(vals)<2:return None
    arr=np.asarray(vals,dtype=float)
    med=np.median(rng.choice(arr,(N,len(arr)),replace=True),axis=1)
    return [float(z) for z in np.quantile(med,[.025,.975])]


def group_stats(records,rng):
    if not records:return {"localities":0}
    keys=("host_report_alien_share","flora_all_alien_share",
          "flora_family_mixture_expected_share","host_minus_flora_all",
          "host_minus_flora_family")
    x={"localities":len(records),"HOSTS_report_rows":sum(r["n"] for r in records),
       "locality_names":[r["location"] for r in records],
       "record_weighted_host_alien_share":sum(r["introduced_host_reports"] for r in records)/sum(r["n"] for r in records)}
    for k in keys:
        v=[r[k] for r in records]
        x["median_"+k]=float(np.median(v))
        x["median_"+k+"_bootstrap95_descriptive"]=quantile_ci(v,rng,N=2999)
    return x


def compare(a,b,key,rng):
    if min(len(a),len(b))<3:return {"status":"NOT_ENOUGH_SOURCE_LOCALITIES","n":[len(a),len(b)]}
    x=[r[key] for r in a];y=[r[key] for r in b]
    vals=[]
    for _ in range(4999):
        p=np.median(rng.choice(x,len(x),replace=True))
        q=np.median(rng.choice(y,len(y),replace=True))
        vals.append(p-q)
    test=mannwhitneyu(x,y,method="auto",alternative="two-sided")
    return {"n":[len(a),len(b)],
            "median_difference":float(np.median(x)-np.median(y)),
            "descriptive_locality_bootstrap95":[float(z) for z in np.quantile(vals,[.025,.975])],
            "exploratory_mannwhitney_p":float(test.pvalue),
            "test_unit":"HOSTS original exact locality; geography and host reporting correlated"}


def report(rows,mode,cohort,threshold):
    rng=np.random.default_rng(20261116)
    a=[r for r in rows if r["is_island"]]
    b=[r for r in rows if not r["is_island"]]
    o=[r for r in a if r["origin_class"]=="oceanic"]
    c=[r for r in a if r["origin_class"]=="continental"]
    res={"mode":mode,"cohort":cohort,"min_reports":threshold,
         "source_localities":len(rows),"source_host_rows":sum(r["n"] for r in rows),
         "island":group_stats(a,rng),"mainland":group_stats(b,rng),
         "oceanic":group_stats(o,rng),"continental":group_stats(c,rng),
         "island_vs_mainland":{
            "raw_report_share":compare(a,b,"host_report_alien_share",rng),
            "whole_flora_adjusted_excess":compare(a,b,"host_minus_flora_all",rng),
            "family_flora_adjusted_excess":compare(a,b,"host_minus_flora_family",rng),
         },
         "oceanic_vs_continental":{
            "family_flora_adjusted_excess":compare(o,c,"host_minus_flora_family",rng),
         },
        }
    res["condition_on_ports_area_reporting_effort"]={
        "family_flora_adjusted":unweighted_lm(
            [{**r,"family_residual":r["host_minus_flora_family"]} for r in rows],
            response="family_residual",include_geology=False,bootstrap=1499)
    }
    return res


def main():
    p=argparse.ArgumentParser()
    for name in ("hosts","host-resolution","wcvp-host-origin","locations","fao",
                 "ports","areas","butterflies239","flora-family","original-island-flora","out-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    src=SimpleNamespace(hosts=a.hosts,host_resolution=a.host_resolution,
           wcvp=a.wcvp_host_origin,locations=a.locations,fao=a.fao,
           ports=a.ports,areas=a.areas,butterflies239=a.butterflies239)
    allrows,locs,coverage,receipt=build_rows(src)
    flor=baseline(a.flora_family)
    control=check_flora(a.flora_family,a.original_island_flora)
    results={}
    records=[]
    skips={}
    for cohort in COHORTS:
        for mode in MODES:
            selected=list(select(allrows,cohort,mode))
            for threshold in (10,20):
                tab,skip=sample(selected,flor,a.host_resolution,threshold,mode)
                key=f"{cohort}__{mode}__min{threshold}"
                results[key]=report(tab,mode,cohort,threshold)
                skips[key]=skip
                if threshold==10 and cohort=="HOSTS_Papilionoidea":
                    records += [{"cropmode":mode,**r} for r in tab]
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"per_locality_family_matched.csv").open("w",newline="",encoding="utf8") as f:
        if not records:raise RuntimeError("No family-matched HOSTS localities")
        w=csv.DictWriter(f,fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)
    audit={
        "status":"CONTEMPORARY_FLORISTIC_BASELINE_NOT_SELECTIVE_FEEDING",
        "estimand":"HOSTS locality-wise introduced host-report fraction minus accepted WCVP regional whole-flora or report-family-matched introduced share",
        "source_base":receipt,"source_coverage":coverage,
        "accepted_WCVP_all_species_source_controls":control,
        "flora_region_family_source":"WCVP accepted plant family-by-region; same FAO exclusion rules; no plant abundances",
        "source_locality_count":len(locs),
        "method":"one exact original location and one WGSRPD3 code, 10+ origin-matched literature rows; all-flora and within-reported-family expected alien fraction",
        "scopes":results,"exclusion_counts":skips,
        "caveats":[
            "Source HOSTS interaction bibliography and original potential resource study reuse the same host catalogue, not independent field observations.",
            "An all-vascular plant flora baseline is NOT the set of compatible larval hosts or its abundance.",
            "Family-matched WCVP accepted flora is not a null for actual feeding preference or realized encounter; availability and phylogenetic relationships within families differ.",
            "Same contemporary WCVP source underlies host origin and flora baselines, making the comparison algebraically and taxonomically dependent.",
            "No field observations of wild larvae, feeding, oviposition, survival, local recruitment or shipping transport were sourced.",
            "Exact single-code locations only; composite/approximate localities, host names with unknown origin or missing WCVP family are excluded.",
            "FAO crop exclusion is taxonomic; not independent evidence of habitat cultivation at a locality.",
            "WPI port counts and island geology are exploratory geographic predictors, not historical host-propagule transport or site-level accessibility.",
            "Source-locality bootstrap and Mann-Whitney p-values are descriptive and do not resolve shared insect, plant or literature publication dependence.",
        ]}
    (a.out_dir/"host_report_minus_family_flora_baseline.json").write_text(
        json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    concise={}
    for key,x in results.items():
        if x["min_reports"]!=10:continue
        concise[key]={
          "n":x["source_localities"],"host_rows":x["source_host_rows"],
          "island":{k:v for k,v in x["island"].items() if k.startswith("median_") or k=="localities"},
          "mainland":{k:v for k,v in x["mainland"].items() if k.startswith("median_") or k=="localities"},
          "contrasts":x["island_vs_mainland"],
          "ports":x["condition_on_ports_area_reporting_effort"]["family_flora_adjusted"],
          "oceanic":x["oceanic"]["localities"],"continental":x["continental"]["localities"],
        }
    print("HOSTS_VS_WCVP_FLORISTIC_BASELINE",json.dumps({
        "source_control":control,"scopes":concise},ensure_ascii=False),flush=True)


if __name__=="__main__":main()

#!/usr/bin/env python3
"""Conservative screen of documented migrants against island GBIF record gaps.

Source Chowdhury et al. 2025 GBIF request contains published-author selected
migrant taxonKeys derived from Chowdhury et al. 2021 review. Exact names only.
Unmatched species are UNKNOWN, NEVER experimentally classified nonmigrants.

This is NOT host-mediated movement or port causality: 2001-2024 pooled flora,
modern global host geography and GBIF first records cannot identify movement.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCE={
 "migration_crosswalk":(
    "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
    "439dd1034ca16f96ff94222695d98dc39747c2e4"
 ),
 "gap_species":(
    "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/matched_species_screen.csv",
    "eb316611c58d5d08e37345f9f944a0dbf170a41e"
 ),
 "fiveyear_events":(
    "provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv",
    "efbc770469a720e2dea82478db2f6196ba03b353"
 )
}
CANDIDATE="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"

def read(name):
    path,expected=SOURCE[name]
    raw=(ROOT/path).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if actual!=expected:
        raise ValueError(f"Published source changed: {path} [{actual}]")
    return list(csv.DictReader(raw.decode("utf-8").splitlines()))

def ranks(a):
    o=sorted(enumerate(a),key=lambda x:x[1])
    ans=[0.0]*len(a);i=0
    while i<len(o):
        j=i+1
        while j<len(o) and o[j][1]==o[i][1]:
            j+=1
        for k in range(i,j):ans[o[k][0]]=(i+j+1)/2
        i=j
    return ans

def corr(a,b):
    if len(a)<3:return None
    avg=sum(a)/len(a);bvg=sum(b)/len(b)
    xx=sum((x-avg)**2 for x in a)
    yy=sum((x-bvg)**2 for x in b)
    return sum((x-avg)*(y-bvg) for x,y in zip(a,b))/math.sqrt(xx*yy) if xx and yy else None

def spearman(a,b):
    return corr(ranks(a),ranks(b))

def quantile(a,q):
    z=sorted(a)
    i=(len(z)-1)*q
    lo=int(i);hi=math.ceil(i)
    return z[lo]+(z[hi]-z[lo])*(i-lo)

def metrics(rows):
    y=[int(x["far500"])/int(x["n"]) for x in rows]
    z=[int(x["far100"])/int(x["n"]) for x in rows]
    v={
        "host_fraction_vs_500km_gap":spearman([float(x["host"]) for x in rows],y),
        "candidate_port_count_vs_500km_gap":spearman([float(x["port"]) for x in rows],y),
        "wing_size_vs_500km_gap":spearman([float(x["wing"]) for x in rows],y),
        "host_fraction_vs_100km_gap":spearman([float(x["host"]) for x in rows],z),
        "candidate_port_count_vs_100km_gap":spearman([float(x["port"]) for x in rows],z),
    }
    return v

def summarise(rows):
    n=sum(int(x["n"]) for x in rows)
    far=sum(int(x["far500"]) for x in rows)
    return {
        "species":len(rows),
        "first_indexed_events_with_before_after_observation_opportunity":n,
        "prior_record_gaps_ge500km":far,
        "event_weighted_gap_share":far/n if n else None,
        "equal_species_mean_gap_share":(
            sum(int(x["far500"])/int(x["n"]) for x in rows)/len(rows) if rows else None),
        "correlations":metrics(rows),
    }

def analyze(draws=1500):
    cw=read("migration_crosswalk")
    gap=read("gap_species")
    ev=read("fiveyear_events")
    if (len(cw),len(gap),len(ev))!=(225,193,5665):
        raise ValueError("Candidate, gap or butterfly-event counts changed")
    migrants={x["species"] for x in cw if x["author_2025_GBIF_requested_taxon_key_proxy"]==CANDIDATE}
    assert len(migrants)==111
    strong=[x for x in gap if x["species"] in migrants]
    unknown=[x for x in gap if x["species"] not in migrants]
    assert len(strong)==106 and len(unknown)==87

    rng=random.Random(20261031)
    tracked=("host_fraction_vs_500km_gap","candidate_port_count_vs_500km_gap","wing_size_vs_500km_gap")
    vectors={v:[] for v in tracked}
    for _ in range(draws):
        sample=[strong[rng.randrange(len(strong))] for _k in range(len(strong))]
        m=metrics(sample)
        for k in tracked:
            if m[k] is not None and math.isfinite(m[k]):
                vectors[k].append(m[k])
    bs={k:{"successful_draws":len(v),
           "species_bootstrap_95":[quantile(v,.025),quantile(v,.975)]}
        for k,v in vectors.items()}

    # Independent source of 2011-2017 butterfly first-index events, but its
    # pooled 2001-2024 botanical host categorisation is *not* historically
    # independent. Never infer survival / establishment from rerecording alone.
    statuses=("ALIEN_CANDIDATE_RECORDED","NATIVE_ONLY_CLASSIFIED","NO_MATCHING_HOST_RECORD")
    sub=[x for x in ev if x["species"] in migrants]
    rates={}
    group=defaultdict(lambda:{"a":[],"n":[]})
    for status in statuses:
        q=[x for x in sub if x["host_group"]==status]
        rates[status]={"events":len(q),"species":len(set(x["species"] for x in q)),
                       "rerecord_within5_years":sum(x["any_followup_record"]=="True" for x in q)}
        rates[status]["rerecord_fraction"]=rates[status]["rerecord_within5_years"]/len(q) if q else None
        if status in statuses[:2]:
            code="a" if status==statuses[0] else "n"
            for x in q:group[x["species"]][code].append(int(x["any_followup_record"]=="True"))
    matched=[]
    for name,v in group.items():
        if v["a"] and v["n"]:
            da=sum(v["a"])/len(v["a"])
            dn=sum(v["n"])/len(v["n"])
            matched.append({"species":name,"n_alien_candidate":len(v["a"]),
                            "n_native_only":len(v["n"]),"difference":da-dn})
    assert len(matched)==12
    rate_equal=sum(x["difference"] for x in matched)/len(matched)

    return {
      "schema":"chocho_independently_documented_migrant_species_host_gap_screen_v1",
      "status":"DESCRIPTIVE_MIGRANT_SUBSET_WEAK_UNSTABLE_HOST_PROXY_ASSOCIATION",
      "source_git_blobs":{k:{"path":path,"sha":sha} for k,(path,sha) in SOURCE.items()},
      "source_migration_taxonomy":"exact canonical names from 2025 author GBIF predicate, derived from 2021 review; indirect, not reviewed S1",
      "migration_classification_warning":"not listed = migration status unknown, NOT nonmigrant",
      "n_source_225":len(cw),
      "n_source_documented_candidate":len(migrants),
      "all_193_screened":len(gap),
      "documented_subset":summarise(strong),
      "non_listed_unknown_coverage_NOT_A_CONTROL":summarise(unknown),
      "min_events_sensitivity":{
        str(n):summarise([x for x in strong if int(x["n"])>=n])
        for n in (5,10,20)
      },
      "bootstrap_seed":20261031,
      "bootstrap_draws":draws,
      "species_resampled_uncertainty":bs,
      "retrospective_5_year_followup":{
         "documented_migrants_all_first_indexed_records":len(sub),
         "pooled_host_status_groups":rates,
         "same_migrant_species_alien_vs_native_comparable_species":len(matched),
         "same_species_equal_mean_reindex_fraction_difference":rate_equal,
         "same_species_direction":{
           "positive":sum(x["difference"]>0 for x in matched),
           "zero":sum(x["difference"]==0 for x in matched),
           "negative":sum(x["difference"]<0 for x in matched)
         }
      },
      "limits":[
         "Independent migratory status is positive-only and indirectly recovered from author GBIF request; not the review's S1.",
         "Migration is not demonstrated local long-distance flight, island arrival, or dispersal ability on the study islands.",
         "The 193-species gap response is distance to earlier GBIF records, not true butterfly travel.",
         "2001-2024 port botanical host records may postdate 2011-2022 butterfly occurrences.",
         "Global host introduced share not observed historically on target islands.",
         "Repeated recording does not demonstrate breeding, establishment, or larval feeding.",
         "Bootstrap resamples species but ignores phylogeny, island/shared-region effects and measurement error.",
         "Not-listed/unknown is NOT a no-migration or no-host control.",
         "Original species×island×year zero-event risk set and year-specific port/nonport plant flora unavailable.",
      ]
    }

def main():
    out=analyze()
    (HERE/"documented_migrant_host_gap_results.json").write_text(
       json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:out[k] for k in ("status","n_source_documented_candidate","documented_subset",
        "min_events_sensitivity","retrospective_5_year_followup")},ensure_ascii=False,indent=2))

if __name__=="__main__":main()

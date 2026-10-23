#!/usr/bin/env python3
"""Exact country-stratified degree-constrained host-IDENTITY reference for GBIF GADM.

Pure DATA-GEOGRAPHY / LABEL-ASSIGNMENT comparator. Source butterfly↔plant
links are real accepted HOSTS/WCVP candidate identities; random alternatives
can be biological NON-hosts. NEVER interpret the conditional tail as a real
ecological p-value, larval feeding or a butterfly population response.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

PLANT_SCHEMA="chocho_all26_gbif_plant_gadm1_feasibility_20261023_v1"
OUTPUT_SCHEMA="chocho_all26_exact_country_degree_host_identity_reference_v1"
CENSUS=(26,24,44,46,19,231,453)
MAX_GRAPHS_PER_COUNTRY=200000
RECORD_GATE=5


def source_sha(raw):
    return hashlib.sha256(raw).hexdigest()


def build_countries(host, butterfly):
    if host.get("schema")!=PLANT_SCHEMA or host.get("status")!="ALL26_BOTANICAL_GADM_QUERIES_COMPLETE":
        raise ValueError("All source botanical country facets must be complete")
    if (host.get("added_butterfly_country_cases"),host.get("taxa_attempted"),
        host.get("plant_country_source_pairs"),
        sum(len(p["linked_butterflies"]) for p in host["plant_country_results"]),
        host.get("cases_with_butterfly_GADM1_overlap_with_indexed_host"))!=(26,24,44,46,19):
        raise ValueError("26-case original source botanical matrix changed")
    if butterfly.get("status")!="COMPLETE_48" or len(butterfly.get("observations",[]))!=48:
        raise ValueError("All 48 original focal butterfly GADM1 cells missing")
    bf={(r["species"],r["country"]):r for r in butterfly["observations"]}
    if len(bf)!=48:
        raise ValueError("Duplicate source butterfly-country")
    cases=defaultdict(list)
    for r in host["butterfly_country_results"]:
        if not r["source_all_plant_facets_complete"]:
            raise ValueError("Unresolved source plant facet not biological negative")
        key=(r["butterfly"],r["country"])
        if key not in bf:
            raise ValueError("Original butterfly GADM1 source case missing")
        gadm=bf[key]["facets"]["gadm1"]["categories"]
        if sum(gadm.values())!=r["all_butterfly_mapped_gadm1_records"]:
            raise ValueError("Focal GADM1 counts drifted")
        cases[r["country"]].append({
            "butterfly":r["butterfly"],
            "candidate_k":r["source_candidate_host_species"],
            "gadm":gadm,"N":sum(gadm.values()),
            "observed_count":r["butterfly_records_in_candidate_plant_positive_gadm1"]
        })
    plants=defaultdict(list)
    for p in host["plant_country_results"]:
        if p["status"]!="COMPLETE_INDEXED_PLANT_FACET":
            raise ValueError("Plant source GADM facet not complete")
        if p["indexed_country_plant_records"] < sum(p["gadm1_counts"].values()):
            raise ValueError("Plant facet total inconsistent")
        plants[p["country"]].append({
            "accepted_id":p["wcvp_accepted_id"],
            "name":p["plant_name"],
            "links":frozenset(p["linked_butterflies"]),
            "gadm":frozenset(k for k,v in p["gadm1_counts"].items() if v>0)
        })
    if set(cases)!=set(plants) or len(cases)!=11:
        raise ValueError("Source host country population missing or extra")
    for ct in cases:
        b=cases[ct]
        p=plants[ct]
        if len({x["butterfly"] for x in b})!=len(b) or len({x["accepted_id"] for x in p})!=len(p):
            raise ValueError("Duplicated country source taxa")
        for row in b:
            obs={x["accepted_id"] for x in p if row["butterfly"] in x["links"]}
            if len(obs)!=row["candidate_k"]:
                raise ValueError("Original butterfly candidate-host row degree changed")
            united=set.union(set(),*(set(x["gadm"]) for x in p if row["butterfly"] in x["links"]))
            check=sum(v for k,v in row["gadm"].items() if k in united)
            if check!=row["observed_count"]:
                raise ValueError("Archived butterfly plant-GADM observation could not be reconstructed")
        for hostp in p:
            if not hostp["links"] <= {x["butterfly"] for x in b}:
                raise ValueError("Plant column host source identity missing")
    if sum(x["N"] for rows in cases.values() for x in rows)!=453:
        raise ValueError("Original 453 focal GADM records changed")
    if sum(x["observed_count"] for rows in cases.values() for x in rows)!=231:
        raise ValueError("Original 231 focal host-country co-detections changed")
    return cases,plants


def enumerate_graphs(butterflies,plants):
    """Uniform enumeration of *simple* bipartite country graphs at fixed degrees.

    Each row chooses a subset of remaining plant columns, hence no two
    representations for the same graph. Exact finite counts, not MCMC.
    """
    butterflies=sorted(butterflies,key=lambda r:r["butterfly"])
    plants=sorted(plants,key=lambda r:r["accepted_id"])
    rows=[r["candidate_k"] for r in butterflies]
    cols=[len(p["links"]) for p in plants]
    if sum(rows)!=sum(cols) or max(cols,default=0)>len(rows) or max(rows,default=0)>len(cols):
        raise ValueError("Bipartite source degrees inconsistent")
    index={p["accepted_id"]:i for i,p in enumerate(plants)}
    obs_sel=[tuple(sorted(index[p["accepted_id"]] for p in plants if b["butterfly"] in p["links"])) for b in butterflies]
    observed_all=observed_ge5=0
    for b,inds in zip(butterflies,obs_sel):
        marked=set().union(*(plants[j]["gadm"] for j in inds))
        s=sum(n for k,n in b["gadm"].items() if k in marked)
        observed_all+=s
        if b["N"]>=RECORD_GATE:observed_ge5+=s
    joint=Counter()
    equal_case_sum=Fraction(0)
    remaining=cols[:]
    current=[]
    evaluated=0
    def walk(i):
        nonlocal evaluated,equal_case_sum
        if i==len(rows):
            if any(remaining):
                return
            evaluated+=1
            if evaluated>MAX_GRAPHS_PER_COUNTRY:
                raise RuntimeError("Too many possible fixed-margin country graphs")
            total=large=0
            case_fraction=Fraction(0)
            for b,subset in zip(butterflies,current):
                together=set().union(*(plants[j]["gadm"] for j in subset))
                n=sum(v for unit,v in b["gadm"].items() if unit in together)
                total+=n
                case_fraction+=Fraction(n,b["N"])
                if b["N"]>=RECORD_GATE:large+=n
            joint[total,large]+=1
            equal_case_sum+=case_fraction
            return
        eligible=[j for j,n in enumerate(remaining) if n>0]
        k=rows[i]
        if len(eligible)<k:
            return
        for choice in combinations(eligible,k):
            for j in choice:remaining[j]-=1
            remaining_rows=len(rows)-i-1
            if sum(remaining)==sum(rows[i+1:]) and all(n<=remaining_rows for n in remaining):
                current.append(choice)
                walk(i+1)
                current.pop()
            for j in choice:remaining[j]+=1
    walk(0)
    if evaluated==0 or sum(joint.values())!=evaluated:
        raise ValueError("No valid country host identity graph")
    if joint[(observed_all,observed_ge5)]<1:
        raise ValueError("Actual plant/butterfly source graph absent from enumerated exact states")
    return {"states":evaluated,"joint":joint,
            "observed_all":observed_all,"observed_n5":observed_ge5,
            "butterflies":len(butterflies),
            "mapped_butterfly_records":sum(r["N"] for r in butterflies),
            "mapped_butterfly_records_n5":sum(r["N"] for r in butterflies if r["N"]>=RECORD_GATE),
            "observed_sum_case_fractions":sum(
                 (Fraction(r["observed_count"],r["N"]) for r in butterflies),Fraction(0)),
            "expected_sum_case_fractions":equal_case_sum/evaluated}


def exact_pmf(country_states,idx):
    dist=Counter({0:1})
    for c in country_states:
        bycountry=Counter()
        for (a,b),n in c["joint"].items():
            bycountry[(a,b)[idx]]+=n
        updated=Counter()
        for s,v in dist.items():
            for score,n in bycountry.items():
                updated[s+score]+=v*n
        dist=updated
    expected=Fraction(sum(k*v for k,v in dist.items()),sum(dist.values()))
    return dist,expected


def analyze(source,butterfly):
    cases,plants=build_countries(source,butterfly)
    countries=[]
    for ct in sorted(cases):
        g=enumerate_graphs(cases[ct],plants[ct])
        marginal=Counter()
        for (score,n5),count in g["joint"].items():
            marginal[score]+=count
        expected=Fraction(sum(k*v for k,v in marginal.items()),g["states"])
        countries.append({**g,"country":ct,"expected_all":expected,
                          "variable_graph_reference":len(marginal)>1,
                          "unique_exact_T_values":len(marginal),
                          "candidate_plant_taxa":len(plants[ct])})
    total_graphs=math.prod(c["states"] for c in countries)
    if all(c["states"]==1 for c in countries):
        raise ValueError("No swappable host-plant identity: null not informative")
    d,E=exact_pmf(countries,0)
    d5,E5=exact_pmf(countries,1)
    if sum(d.values())!=total_graphs or sum(d5.values())!=total_graphs:
        raise ValueError("Exact country-label distribution lost mass")
    O=sum(c["observed_all"] for c in countries)
    O5=sum(c["observed_n5"] for c in countries)
    N=sum(c["mapped_butterfly_records"] for c in countries)
    N5=sum(c["mapped_butterfly_records_n5"] for c in countries)
    if (O,N)!=(231,453):
        raise ValueError("Original source overlap/nominator changed")
    obs_balanced=sum((c["observed_sum_case_fractions"] for c in countries),Fraction(0))/26
    expected_balanced=sum((c["expected_sum_case_fractions"] for c in countries),Fraction(0))/26
    tail=Fraction(sum(v for k,v in d.items() if k>=O),total_graphs)
    tail5=Fraction(sum(v for k,v in d5.items() if k>=O5),total_graphs)
    leave=[]
    for c in countries:
        this_n=N-c["mapped_butterfly_records"]
        this_o=O-c["observed_all"]
        this_e=E-c["expected_all"]
        leave.append({
            "omitted_country":c["country"],"remaining_butterfly_records":this_n,
            "remaining_observed_host_GADM_records":this_o,
            "remaining_expected_records":round(float(this_e),8),
            "observed_minus_reference_pp":round(100*float((this_o-this_e)/this_n),6)
        })
    return {
        "schema":OUTPUT_SCHEMA,"status":"COMPLETE_EXACT_POSTHOC_HOST_IDENTITY_REFERENCE_NOT_BIOLOGICAL_TEST",
        "source_countries":len(countries),"source_butterfly_country_cases":26,
        "source_butterfly_GADM_records":N,"source_host_positive_GADM_butterfly_records":O,
        "source_host_positive_GADM_focal_fraction":round(O/N,9),
        "exact_expected_host_positive_GADM_records":round(float(E),10),
        "observed_minus_identity_reference_pp":round(100*float((O-E)/N),7),
        "reference_upper_tail_NOT_BIOLOGICAL_P":round(float(tail),10),
        "total_possible_country_degree_fixed_graphs":str(total_graphs),
        "countries_with_more_than_one_graph":sum(c["states"]>1 for c in countries),
        "countries_with_nonconstant_T":sum(c["variable_graph_reference"] for c in countries),
        "source_butterfly_records_min5":{"cases":sum(r["N"]>=RECORD_GATE for group in cases.values() for r in group),
           "record_min":5,"N":N5,"O":O5,"exact_expected":round(float(E5),10),
           "excess_pp":round(100*float((O5-E5)/N5),7),
           "reference_upper_tail_NOT_BIOLOGICAL_P":round(float(tail5),10)},
        "balanced_case_fraction":{"observed":round(float(obs_balanced),9),
              "reference_expected":round(float(expected_balanced),9),
              "observed_minus_reference_pp":round(100*float(obs_balanced-expected_balanced),7)},
        "country_source_statistics":[{
             "country":c["country"],"butterflies":c["butterflies"],
             "candidate_plant_taxa":c["candidate_plant_taxa"],
             "degree_fixed_graph_states":c["states"],
             "geography_changes_under_graph_permutation":c["variable_graph_reference"],
             "country_N":c["mapped_butterfly_records"],
             "country_observed":c["observed_all"],
             "country_expected":round(float(c["expected_all"]),10),
             "country_observed_minus_expected":round(float(c["observed_all"]-c["expected_all"]),10),
             } for c in countries],
        "leave_one_country_out_posthoc":leave,
        "exact_country_convolution_T_counts":{str(k):str(v) for k,v in sorted(d.items())},
        "exact_country_convolution_T_min5_counts":{str(k):str(v) for k,v in sorted(d5.items())},
        "caution":"Plant swaps may link a butterfly to a botanical NONHOST. Finite label-allocation tails are not ecological probabilities, reproductive success tests, confidence levels or evidence of host-mediated population growth.",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-botanical",required=True,type=Path)
    p.add_argument("--original-butterfly-48",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    raw=a.source_botanical.read_bytes()
    butterfly_raw=a.original_butterfly_48.read_bytes()
    output=analyze(json.loads(raw),json.loads(butterfly_raw))
    output["source_SHA256"]={"botanical_gadm1":source_sha(raw),
                            "butterfly_original48":source_sha(butterfly_raw)}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("ALL26_HOST_IDENTITY_REFERENCE="+json.dumps(
        {k:v for k,v in output.items() if k not in
         ("exact_country_convolution_T_counts","exact_country_convolution_T_min5_counts")},
        separators=(",",":")),flush=True)


if __name__=="__main__":
    main()

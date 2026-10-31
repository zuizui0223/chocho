#!/usr/bin/env python3
"""Migrant-subset island host similarity and prevalence-preserving label null.

This analysis is retrospective. Island-host Jaccard uses plants recorded in
2001-2024, potentially AFTER the butterfly outcomes (2011-2017). Positive
author-selected migration status derives indirectly from the 2025 GBIF
request based on Chowdhury et al., not individual butterflies flying to
those islands; not-listed status is UNKNOWN, not nonmigrant.

Permutation preserves only each species' number of source first-indexed
islands within broad prevalence strata. It does not establish geographic,
taxonomic, phylogenetic, observer-effort or shipping exchangeability.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
FILES={
    "migration":(
        "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
        "439dd1034ca16f96ff94222695d98dc39747c2e4"),
    "islands":(
        "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv",
        "18cecde5bf3066f276380167f451f1fa5cc1caa1"),
    "pairs":(
        "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv",
        "5663d2ffe42419f5c4496b92c9f3bbf76175907b"),
    "events":(
        "provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv",
        "efbc770469a720e2dea82478db2f6196ba03b353")
}
KNOWN="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"
UNKNOWN="NOT_LISTED_OR_TAXON_UNRESOLVED_UNKNOWN"
BOOT_SEED=20261031


def read(name):
    path,expected=FILES[name]
    raw=(ROOT/path).read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if blob!=expected:
        raise ValueError("Source Git blob changed: "+path)
    return list(csv.DictReader(io.StringIO(raw.decode("utf8"))))


def midrank(arr):
    pairs=sorted(enumerate(arr),key=lambda item:item[1])
    out=[0.0]*len(pairs)
    i=0
    while i<len(pairs):
        j=i+1
        while j<len(pairs) and pairs[i][1]==pairs[j][1]:
            j+=1
        for k in range(i,j):
            out[pairs[k][0]]=(i+j+1)/2
        i=j
    return out


def pearson(a,b):
    if len(a)<3:
        return None
    ma=sum(a)/len(a); mb=sum(b)/len(b)
    xy=sum((x-ma)*(y-mb) for x,y in zip(a,b))
    xx=sum((x-ma)**2 for x in a)
    yy=sum((y-mb)**2 for y in b)
    return xy/math.sqrt(xx*yy) if xx>0 and yy>0 else None


def spearman(a,b):
    return pearson(midrank(a),midrank(b))


def solve(A,b):
    n=len(b)
    x=[A[i][:]+[b[i]] for i in range(n)]
    for k in range(n):
        best=max(range(k,n),key=lambda i:abs(x[i][k]))
        if abs(x[best][k])<1e-12:
            return None
        x[k],x[best]=x[best],x[k]
        pivot=x[k][k]
        for j in range(k,n+1):
            x[k][j]/=pivot
        for i in range(n):
            if i==k:continue
            v=x[i][k]
            for j in range(k,n+1):
                x[i][j]-=v*x[k][j]
    return [r[n] for r in x]


def residual(y,X):
    k=len(X[0])
    AA=[[sum(row[i]*row[j] for row in X) for j in range(k)] for i in range(k)]
    bb=[sum(row[i]*v for row,v in zip(X,y)) for i in range(k)]
    beta=solve(AA,bb)
    return ([v-sum(z*x for z,x in zip(row,beta)) for row,v in zip(X,y)]
            if beta is not None else None)


def partial_rank(host,obs,distance,prior):
    h,y,d,p=map(midrank,(host,obs,distance,prior))
    X=[[1.0,d[i],p[i]] for i in range(len(h))]
    x=residual(h,X)
    z=residual(y,X)
    return pearson(x,z) if x is not None and z is not None else None


def jaccard(a,b):
    union=a|b
    if not union:
        return None
    return len(a&b)/len(union)


def base_tables():
    islands=[r["landmass_id"] for r in read("islands") if r["eligible"]=="True"]
    if len(islands)!=51 or len(set(islands))!=51:
        raise ValueError("Original 51 island support changed")
    idset=set(islands)
    pairs=[]
    for r in read("pairs"):
        a,b=r["island_a"],r["island_b"]
        if a not in idset or b not in idset or a==b:
            raise ValueError("Island pair outside original eligible islands")
        pairs.append({
            "a":a,"b":b,
            "host":float(r["host_jaccard"]),
            "prior":float(r["earlier_butterfly_jaccard"]),
            "distance":float(r["port_distance_km"]),
            "all_new_2011_2022":float(r["new_butterfly_jaccard"])
        })
    assert len(pairs)==1275 and len({(r["a"],r["b"]) for r in pairs})==1275
    migration={}
    for r in read("migration"):
        x=r["author_2025_GBIF_requested_taxon_key_proxy"]
        if x not in (KNOWN,UNKNOWN):
            raise ValueError("Unknown roster category")
        key=r["speciesKey"]
        if key in migration:raise ValueError("Duplicate GBIF speciesKey in roster")
        migration[key]=x==KNOWN
    assert len(migration)==225
    first=[r for r in read("events") if r["landmass_id"] in idset]
    species_islands={}
    for r in first:
        key=r["speciesKey"]
        if key not in migration:
            raise ValueError("Unknown butterfly species key")
        if key not in species_islands:
            species_islands[key]=set()
        if r["landmass_id"] in species_islands[key]:
            raise ValueError("Repeated butterfly first island index")
        species_islands[key].add(r["landmass_id"])
    assert len(first)==302 and len(species_islands)==106
    marked={s for s in species_islands if migration[s]}
    assert len(marked)==63
    return islands,pairs,species_islands,marked


def score(ids,pairs,species_islands,marked):
    by_island={i:set() for i in ids}
    for sp in marked:
        for island in species_islands[sp]:
            if island in by_island:
                by_island[island].add(sp)
    filtered=[]
    for p in pairs:
        value=jaccard(by_island[p["a"]],by_island[p["b"]])
        if value is not None:
            filtered.append((p,value))
    hosts=[p["host"] for p,_ in filtered]
    obs=[v for _,v in filtered]
    dist=[p["distance"] for p,_ in filtered]
    prior=[p["prior"] for p,_ in filtered]
    return {
        "n_pairs":len(filtered),
        "jaccard_zero_pairs":sum(v==0 for v in obs),
        "spearman_host_vs_butterfly":spearman(hosts,obs),
        "partial_rank_distance_prior_butterfly":partial_rank(hosts,obs,dist,prior)
    }


def binning(islands_by_species):
    buckets={}
    for species,where in islands_by_species.items():
        n=len(where)
        key="1" if n==1 else ("2" if n==2 else ("3-4" if n<=4 else "5plus"))
        buckets.setdefault(key,[]).append(species)
    return buckets


def rng_choice(state,n):
    state=(1664525*state+1013904223)&0xffffffff
    return state,state%n


def qtls(arr):
    arr=sorted(arr)
    def q(p):
        k=(len(arr)-1)*p
        a,b=int(math.floor(k)),int(math.ceil(k))
        return arr[a]+(arr[b]-arr[a])*(k-a)
    return [q(.025),q(.975)]


def analyze(draws=999):
    islands,pairs,species_islands,observed_marked=base_tables()
    pool=set(species_islands)
    assert len(pool)==106
    score_m=score(islands,pairs,species_islands,observed_marked)
    score_u=score(islands,pairs,species_islands,pool-observed_marked)
    score_all=score(islands,pairs,species_islands,pool)
    baseline={
        "n_pairs":len(pairs),
        "spearman_host_vs_butterfly":spearman(
            [p["host"] for p in pairs],
            [p["all_new_2011_2022"] for p in pairs]),
        "partial_rank_distance_prior_butterfly":partial_rank(
            [p["host"] for p in pairs],
            [p["all_new_2011_2022"] for p in pairs],
            [p["distance"] for p in pairs],
            [p["prior"] for p in pairs])
    }
    assert abs(baseline["spearman_host_vs_butterfly"]-0.29980576483377025)<1e-10
    assert abs(score_m["spearman_host_vs_butterfly"]-0.2204520274075129)<1e-10

    # Leave-one-island sensitivity, NOT an independence-based confidence interval.
    loio=[
        score([i for i in islands if i!=one],
              [p for p in pairs if one not in (p["a"],p["b"])],
              species_islands,observed_marked)
        for one in islands
    ]
    buckets=binning(species_islands)
    labels={
        name:{"total":len(v),
              "published_migrant_selected":sum(s in observed_marked for s in v)}
        for name,v in buckets.items()
    }
    state=BOOT_SEED
    null=[]
    for _ in range(draws):
        selected=set()
        for name,group in buckets.items():
            v=group[:]
            for j in range(len(v)-1,0,-1):
                state,k=rng_choice(state,j+1)
                v[j],v[k]=v[k],v[j]
            selected.update(v[:labels[name]["published_migrant_selected"]])
        null.append(score(islands,pairs,species_islands,selected))

    nullstats={}
    for metric in ("spearman_host_vs_butterfly","partial_rank_distance_prior_butterfly","n_pairs"):
        values=[x[metric] for x in null]
        nullstats[metric]={
            "mean":sum(values)/len(values),
            "reference95":qtls(values),
            "draws_ge_observed":sum(x>=score_m[metric] for x in values),
            "draws":draws,
            "range":[min(values),max(values)]
        }
    report={
        "schema":"chocho_author_migrant_island_2011_2017_host_network_prevalence_null_v1",
        "status":"POSITIVE_ISLAND_SIMILARITY_NOT_MIGRANT_SPECIFIC_ABOVE_PREVALENCE_NULL",
        "source_git_blobs":{k:{"path":v[0],"sha":v[1]} for k,v in FILES.items()},
        "source_provenance":"2025 author GBIF predicate based on published migratory butterfly review; 2021 Table S1 not obtained",
        "all_225_candidate_species_retained":225,
        "islands":len(islands),
        "island_pairs":len(pairs),
        "butterfly_first_index_2011_2017_island_events":sum(len(v) for v in species_islands.values()),
        "butterfly_species_first_indexed_2011_2017":len(species_islands),
        "author_roster_migrant_selected_species_with_new_island_records":len(observed_marked),
        "not_listed_unknown_species_with_new_records":len(pool-observed_marked),
        "original_all_225_new_2011_2022_baseline":baseline,
        "first_2011_2017_all_species":score_all,
        "first_2011_2017_author_migrant_roster":score_m,
        "first_2011_2017_not_listed_UNKNOWN_not_nonmigrants":score_u,
        "migrant_leave_one_island_sensitivity":{
            "simple_rho_range":[min(q["spearman_host_vs_butterfly"] for q in loio),
                                max(q["spearman_host_vs_butterfly"] for q in loio)],
            "partial_rank_range":[min(q["partial_rank_distance_prior_butterfly"] for q in loio),
                                  max(q["partial_rank_distance_prior_butterfly"] for q in loio)]
        },
        "source_prevalence_bins":labels,
        "permutation":{
            "method":"permute author migration labels within number-of-first-indexed-islands bins (1, 2, 3-4, >=5), keeping bins counts; NOT a geographic/effort causal null",
            "seed":BOOT_SEED,
            "draws":draws,
            "statistics":nullstats,
            "posthoc_descriptive_only":True
        },
        "explicit_limits":[
            "Species migration status in 2025 GBIF predicate is indirect, not original 2021 Table S1.",
            "Nonlisted means unknown, not sedentary/nonmigrant.",
            "Migration in a species somewhere does not establish flight or arrival at these islands.",
            "Butterfly first GBIF index is not first true colonization nor breeding.",
            "Port host composition pools 2001-2024 botanical records, potentially postdating 2011-2017 butterfly indices.",
            "The 51 islands and 1275 pairs are not independent global island samples.",
            "Partial ranks adjust distance and earlier butterfly similarity only; not climate, sampling, phylogeny or port presence.",
            "Migration-label shuffles within prevalence bins do not preserve geography, taxon-family, host geography or source-specific effort.",
            "The 2021-22 new butterfly base and 2011-17 new subset are different time windows.",
            "No evidence of shipborne host introduction, actual larval resource use, establishment or stepping-stone movement."
        ]
    }
    return report


def main():
    result=analyze()
    (HERE/"migrant_island_host_prevalence_null_results.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"n_islands":result["islands"],
                      "observed":result["first_2011_2017_author_migrant_roster"],
                      "permutation":result["permutation"]},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

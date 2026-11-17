#!/usr/bin/env python3
"""Source-locked 225x51 GBIF recording effort / butterfly origin claim gate.

No GBIF non-record is coded as absence or new colonization. Butterfly native
range and endemicity require independent exact-island evidence; WCVP plant
origins and HOSTS native-host geography are NOT butterfly range evidence.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
PIN = "b6a3c98eb76864052d45d520a55b7d830816126c"
DEPTHS = (20, 50, 100, 200)
EVIDENCE_FIELDS = ("speciesKey","landmass_id","origin_status","establishment_status",
                   "endemicity_status","geographic_resolution","source_url",
                   "global_native_range_reviewed")
ALLOWED = ({"native","introduced","vagrant","unknown"},
           {"established","reproducing","casual","unknown"},
           {"single_island_endemic","archipelago_endemic","not_endemic","unassessed"})

def git_blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def load_source(path):
    blob=git_blob(path)
    if blob!=PIN: raise ValueError("Pinned GBIF source SHA drift: "+blob)
    with path.open(newline="",encoding="utf-8-sig") as f:
        rd=csv.DictReader(f)
        must={"species","speciesKey","landmass_id","pre2001_2010_indexed_records",
              "future2011_2017_indexed_records","future2018_2022_indexed_records",
              "first_indexed_in_2011_2017","first_indexed_in_2018_2022"}
        if not must.issubset(rd.fieldnames or []): raise ValueError("Source schema drift")
        rows=[]
        for row in rd:
            a,b,c=[int(row[k]) for k in ("pre2001_2010_indexed_records",
                         "future2011_2017_indexed_records","future2018_2022_indexed_records")]
            if min(a,b,c)<0: raise ValueError("Negative indexed records")
            if int(row["first_indexed_in_2011_2017"])!=int(a==0 and b>0):
                raise ValueError("Early later first-index event inconsistent")
            if int(row["first_indexed_in_2018_2022"])!=int(a==0 and b==0 and c>0):
                raise ValueError("Late later first-index event inconsistent")
            rows.append((row["species"].strip(),int(row["speciesKey"]),int(row["landmass_id"]),a,b,c))
    sp=sorted({r[1] for r in rows}); islands=sorted({r[2] for r in rows})
    if (len(sp),len(islands),len(rows))!=(225,51,11475): raise ValueError("Panel extent drift")
    if len({(r[1],r[2]) for r in rows})!=len(rows): raise ValueError("Duplicated species-island key")
    return rows,sp,islands,blob

def evidence_rows(path,valid):
    if path is None: return {}
    out={}
    with path.open(newline="",encoding="utf-8-sig") as f:
        rd=csv.DictReader(f)
        if not set(EVIDENCE_FIELDS).issubset(rd.fieldnames or []):
            raise ValueError("Curated evidence missing required fields")
        for x in rd:
            key=(int(x["speciesKey"]),int(x["landmass_id"]))
            if key not in valid or key in out: raise ValueError("Unknown or duplicate evidence key")
            o=x["origin_status"].strip(); e=x["establishment_status"].strip()
            n=x["endemicity_status"].strip()
            if o not in ALLOWED[0] or e not in ALLOWED[1] or n not in ALLOWED[2]:
                raise ValueError("Unsupported origin/establishment/endemicity label")
            if o!="unknown" or e!="unknown" or n!="unassessed":
                if x["geographic_resolution"].strip()!="exact_island_landmass" or not x["source_url"].strip().startswith(("https://","http://")):
                    raise ValueError("Evidence must have exact-island geographic scope and source URL")
            if n in {"single_island_endemic","archipelago_endemic","not_endemic"}:
                if x["global_native_range_reviewed"].strip()!="1":
                    raise ValueError("Endemicity requires reviewed complete native range")
            if n in {"single_island_endemic","archipelago_endemic"} and o!="native":
                raise ValueError("Endemicity requires confirmed native origin")
            out[key]={"origin":o,"establishment":e,"endemicity":n,"source_url":x["source_url"].strip()}
    return out

def raw_jaccard(islands,pre,cumulative):
    old=new=0.
    n=up=0
    for pos,i in enumerate(islands):
        for j in islands[pos+1:]:
            a=pre[i];b=pre[j];c=cumulative[i];d=cumulative[j]
            if not a|b or not c|d: continue
            x=len(a&b)/len(a|b);y=len(c&d)/len(c|d)
            old+=x;new+=y;n+=1;up+=y>x
    return {"pairs":n,"early":old/n,"cumulative":new/n,
            "delta":(new-old)/n,"pairs_increased":up}

def include_probability(c,total,m):
    if not c: return 0.
    if total-c<m: return 1.
    miss=1.
    for j in range(m): miss*=(total-c-j)/(total-j)
    return 1.-miss

def weighted_jaccard(p,q):
    both=union=0.
    for x,y in zip(p,q):
        z=x*y;both+=z;union+=x+y-z
    return both/union if union else None

def standardized(islands,species,records,depth):
    totals={i:sum(v[0] for v in records[i].values()) for i in islands}
    pre={};cumulative={}
    for i in islands:
        early=[records[i][s][0] for s in species]
        later=[sum(records[i][s]) for s in species]
        T=totals[i];F=sum(later)
        pre[i]=[include_probability(v,T,depth) for v in early]
        cumulative[i]=[include_probability(v,F,depth) for v in later]
    e=l=0.;n=0
    for pos,i in enumerate(islands):
        for j in islands[pos+1:]:
            x=weighted_jaccard(pre[i],pre[j])
            y=weighted_jaccard(cumulative[i],cumulative[j])
            if x is None or y is None: continue
            e+=x;l+=y;n+=1
    return {"records_per_landmass":depth,"islands":len(islands),"pairs":n,
            "early":e/n,"cumulative":l/n,"delta":(l-e)/n,
            "estimand":"mean of expected_overlap / expected_union across island pairs, not mean E[Jaccard]"}

def analyze(rows,species,islands,blob,evidence):
    records={i:{} for i in islands};pre=defaultdict(set);cum=defaultdict(set)
    early_counts=Counter();new=[]
    for name,s,i,a,b,c in rows:
        records[i][s]=(a,b,c);early_counts[i]+=a
        if a: pre[i].add(s)
        if a+b+c: cum[i].add(s)
        if a==0 and b+c>0:
            status=evidence.get((s,i),{"origin":"unknown","establishment":"unknown",
                                        "endemicity":"unassessed"})
            new.append({"species":name,"speciesKey":s,"landmass_id":i,
                        "first_indexed_window":"2011-2017" if b else "2018-2022",
                        "GBIF_2001_2010_prior_other_butterfly_records":0,
                        "butterfly_origin":status["origin"],
                        "establishment":status["establishment"],
                        "endemicity":status["endemicity"]})
    for k in new:k["GBIF_2001_2010_prior_other_butterfly_records"]=early_counts[k["landmass_id"]]
    dominant=max(islands,key=lambda i:early_counts[i])
    if not all(len(records[i])==len(species) for i in islands): raise ValueError("Unbalanced source")
    depths=[]
    for d in DEPTHS:
        I=[i for i in islands if early_counts[i]>=d]
        depths.append({"min_prior_GBIF_butterfly_records":d,"unthinned":raw_jaccard(I,pre,cum),
                       "equal_record_depth":standardized(I,species,records,d)})
    total=sum(early_counts.values())
    known=sum(x["origin"]!="unknown" for x in evidence.values())
    return {
        "schema":"chocho_butterfly_GBIF_effort_native_endemism_source_gate_20261117_v1",
        "status":"DESCRIPTIVE_RECORD_EFFORT_WITH_INDEPENDENT_BUTTERFLY_ORIGIN_GATE",
        "source_path":REL,"source_git_blob":blob,
        "cohort":{"species":len(species),"landmass_groups":len(islands),
                  "species_island_cells":len(rows),"early_indexed_cells":sum(map(len,pre.values())),
                  "cumulative_indexed_cells":sum(map(len,cum.values())),
                  "later_first_indexed_cells":len(new),"early_total_indexed_records":total,
                  "dominant_landmass_id":dominant,
                  "dominant_landmass_early_indexed_records":early_counts[dominant],
                  "dominant_landmass_record_fraction":early_counts[dominant]/total,
                  "n_early_lt20_landmasses":sum(early_counts[i]<20 for i in islands),
                  "n_first_indices_on_early_lt20_landmasses":sum(early_counts[x["landmass_id"]]<20 for x in new)},
        "unthinned_all51":raw_jaccard(islands,pre,cum),
        "unthinned_without_dominant":raw_jaccard([i for i in islands if i!=dominant],pre,cum),
        "effort_standardized_sensitivity":depths,
        "native_endemism":{"exact_island_curated_evidence_rows":len(evidence),
                           "cells_with_no_independently_sourced_native_label":len(rows)-known,
                           "first_indices_without_origin_label":sum(x["butterfly_origin"]=="unknown" for x in new),
                           "confirmed_native_established_first_indices":sum(x["butterfly_origin"]=="native" and x["establishment"]=="established" for x in new),
                           "verified_single_island_endemic_first_indices":sum(x["endemicity"]=="single_island_endemic" for x in new)},
        "interpretation":["GBIF first index != biological arrival or establishment",
                          "GBIF zero record != verified ecological absence",
                          "Recording counts are imperfect sampling-effort proxies, not survey hours",
                          "HOSTS native host geographic footprint and WCVP plant origin are NOT butterfly native ranges",
                          "Butterfly endemism cannot be identified by GBIF country count or lack of records elsewhere",
                          "A genuine endemic label needs independently reviewed global native range and island scope",
                          "Record-level rarefaction standardizes numbers of GBIF records, NOT observer selection, protocols, or habitat coverage"]
    },new

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--origin-evidence",type=Path)
    p.add_argument("--out-dir",type=Path,default=Path("_butterfly_effort_origin_gate"))
    a=p.parse_args()
    rows,species,islands,blob=load_source(ROOT/REL)
    valid={(r[1],r[2]) for r in rows}
    ev=evidence_rows(a.origin_evidence,valid)
    out,new=analyze(rows,species,islands,blob,ev)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"BUTTERFLY_EFFORT_ORIGIN_GATE_20261117.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
    with (a.out_dir/"LATER_FIRST_INDICES_TO_REVIEW.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(new[0]))
        w.writeheader();w.writerows(sorted(new,key=lambda x:(x["GBIF_2001_2010_prior_other_butterfly_records"],x["landmass_id"],x["speciesKey"])))
    print(json.dumps({"cohort":out["cohort"],"depths":out["effort_standardized_sensitivity"],
                      "origin":out["native_endemism"]},indent=2))

if __name__=="__main__": main()

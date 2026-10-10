#!/usr/bin/env python3
"""Audit frozen delayed-butterfly / shared parasitoid time series support.

This script ONLY inspects source quality, denominator consistency and independent
site-week transitions. It does not fit or select statistical models.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,collections
from pathlib import Path
from urllib.request import Request,urlopen

URL="https://ndownloader.figshare.com/files/29157873"
MD5="5e861203477ce8d2ae3efdfacc7436fc"
META_URL="https://ndownloader.figshare.com/files/29157876"
META_MD5="107ce0df69d9157adabf9572fd94da92"
RESIDENT=["Aglais urticae","Aglais io"]
NEWCOMER="Araschnia levana"
SPECIES_CODES={"au":"Aglais urticae","aio":"Aglais io","alev":"Araschnia levana","va":"Vanessa atalanta"}
def normalized_species(v):
    value=str(v or "").strip().lower()
    if value not in SPECIES_CODES:
        raise ValueError("unrecognized original butterfly taxon code: "+value)
    return SPECIES_CODES[value]
REQUIRED={"year","week_ISO","BMS_id","butterfly_species","nb_larvae_in_lab",
          "larvae_Sturmia_bella","pupae_Sturmia_bella",
          "county","presence_alev","lat4326","instar_rond","total_larvae"}

def retrieve(url,md5):
    with urlopen(Request(url,headers={"User-Agent":"chocho-academic-source-integrity/1.0"}),timeout=35) as r:
        raw=r.read(2_000_000)
    if hashlib.md5(raw).hexdigest()!=md5:
        raise ValueError("original source MD5 mismatch")
    return raw

def integer(s):
    v=str(s or "").strip()
    if v.lower() in ("na","nan","null","","none","-"):
        return None
    f=float(v)
    if not f.is_integer():
        raise ValueError(f"noninteger count: {v}")
    return int(f)

def open_source(raw):
    rr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter=";")
    if not REQUIRED.issubset(rr.fieldnames or []):
        raise ValueError("missing source fields "+repr(sorted(REQUIRED-set(rr.fieldnames or []))))
    rows=list(rr)
    if len(rows)!=1080:
        raise ValueError("not 1080 deposited batch monitoring rows")
    return rows

def source_report(rows):
    species=collections.Counter()
    missing=collections.Counter()
    timekeys=collections.defaultdict(list)
    counters=collections.Counter()
    names=set()
    sites_per_year=collections.defaultdict(set)
    site_counties=collections.defaultdict(set)
    for r in rows:
        original_code=r["butterfly_species"].strip().lower()
        sp=normalized_species(original_code)
        species[original_code]+=1
        for k in ("year","week_ISO","nb_larvae_in_lab","larvae_Sturmia_bella",
                  "pupae_Sturmia_bella","presence_alev","total_larvae","lat4326","instar_rond"):
            if not str(r.get(k) or "").strip() or str(r[k]).strip().lower() in ("na","nan"):
                missing[k]+=1
        id_=r["BMS_id"].strip()
        year=integer(r["year"])
        week=integer(r["week_ISO"])
        n=integer(r["nb_larvae_in_lab"])
        lar=integer(r["larvae_Sturmia_bella"])
        pup=integer(r["pupae_Sturmia_bella"])
        if not id_ or year is None or week is None or not 1<=week<=53:
            counters["invalid_site_week"]+=1
            continue
        timekeys[(id_,year,week)].append((sp,n,lar,pup))
        names.add(id_)
        sites_per_year[year].add(id_)
        site_counties[id_].add(r["county"].strip())
        if sp in RESIDENT:
            counters["resident_batches"]+=1
            if n is None or n<=0 or lar is None or pup is None or lar<0 or pup<0 or lar+pup>n:
                counters["resident_invalid_outcome"]+=1
            else:
                counters["resident_larvae"]+=n
                counters["resident_sturmia"]+=lar+pup
        if sp==NEWCOMER and n is not None and n>0:
            counters["newcomer_batches"]+=1
            counters["newcomer_larvae"]+=n
    events=collections.defaultdict(lambda:{"n_batches":0,"n_larvae":0,"n_sturmia":0,"sites":set()})
    observed_siteweeks=collections.defaultdict(set)
    aleve=collections.Counter()
    for (site,year,week),vals in timekeys.items():
        observed_siteweeks[(site,year)].add(week)
        aleve[(site,year,week)]=any(sp==NEWCOMER and n is not None and n>0 for sp,n,_,_ in vals)
    visited=collections.Counter()
    for (site,year,week),vals in timekeys.items():
        earlier=[w for w in observed_siteweeks[(site,year)] if w<week]
        if not earlier:
            continue
        prev=max(earlier)
        gap=week-prev
        if gap not in (1,2,3):
            visited["disallowed_gap"]+=1
            continue
        p=int(aleve[(site,year,prev)])
        c=int(aleve[(site,year,week)])
        for sp,n,lar,pup in vals:
            if sp not in RESIDENT or n is None or n<=0 or lar is None or pup is None or lar+pup>n:
                continue
            cat=f"previous_{p}_current_{c}"
            a=events[cat]
            a["n_batches"]+=1
            a["n_larvae"]+=n
            a["n_sturmia"]+=lar+pup
            a["sites"].add(site)
            visited["eligible_batches"]+=1
    results={k:{**{kk:vv for kk,vv in v.items() if kk!="sites"},"n_sites":len(v["sites"])} for k,v in events.items()}
    previous_1=[z for k,z in results.items() if k.startswith("previous_1")]
    previous_0=[z for k,z in results.items() if k.startswith("previous_0")]
    prev1=sum(z["n_batches"] for z in previous_1)
    prev0=sum(z["n_batches"] for z in previous_0)
    inf1=sum(z["n_sturmia"] for z in previous_1)
    inf0=sum(z["n_sturmia"] for z in previous_0)
    gates={
        "at_least_10_study_sites":len(names)>=10,
        "at_least_30_prior_presence_batches":prev1>=30,
        "at_least_30_prior_absence_batches":prev0>=30,
        "at_least_10_sturmia_in_both_prior_groups":min(inf1,inf0)>=10,
        "at_least_15_prev_present_now_absent":results.get("previous_1_current_0",{}).get("n_batches",0)>=15,
        "at_least_15_prev_absent_now_present":results.get("previous_0_current_1",{}).get("n_batches",0)>=15,
        "no_invalid_resident_outcomes":counters["resident_invalid_outcome"]==0,
    }
    return {"schema":"chocho_nettle_delayed_enemy_source_support_v01",
            "batch_rows":len(rows),"butterfly_species_labels":dict(species),
            "sites":len(names),
            "original_BMS_id_values":sorted(names),
            "original_BMS_id_count_by_year":{str(k):len(v) for k,v in sites_per_year.items()},
            "county_codes_by_original_BMS_id":{k:sorted(v) for k,v in sorted(site_counties.items())},
            "published_paper_site_count":19,
            "source_code_site_count_mismatch":len(names)!=19,
            "unique_site_year_week":len(timekeys),
            "missing_value_counts":dict(missing),
            "count_integrity":dict(counters),
            "site_week_gaps":dict(visited),
            "previous_by_current_support":results,
            "checks":gates,
            "data_sufficiency_gate_pass":all(gates.values()),
            "effect_estimated":False,
            "not_a_botanical_invasion_test":True,
            "GEB_PR38_unchanged":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    a=p.parse_args()
    original=retrieve(URL,MD5)
    records=open_source(original)
    result=source_report(records)
    result["source_original_md5"]=MD5
    result["source_original_sha256"]=hashlib.sha256(original).hexdigest()
    original_meta=retrieve(META_URL,META_MD5)
    result["source_metadata_sha256"]=hashlib.sha256(original_meta).hexdigest()
    result["metadata_excerpt"]=original_meta.decode("utf-8-sig",errors="replace").splitlines()[:22]
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()

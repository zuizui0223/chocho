#!/usr/bin/env python3
"""Link UK LTMN source XLSX transects to exact original visit histories.

No easting/northing/gridref, transect IDs, original site names or site-year
individual records serialized. Habitat class counts are descriptive, not
observed larval host or proof of human-maintained plants.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math
from pathlib import Path
from urllib.request import Request,urlopen
from audit_uk_ltmn_site_crosswalk_20261028 import URL,SHA,original_rows
from audit_uk_ltmn_historical_event_join_20261028 import get_source,excel_rows,FIELDS_VISITS,parse_event,key_site

def read_crosswalk():
    with urlopen(Request(URL,headers={"User-Agent":"chocho-UK-historical-transect-crosswalk/1.0"}),timeout=45) as res:
        body=res.read(400000)
    if hashlib.sha256(body).hexdigest()!=SHA:raise ValueError("Original LTMN crosswalk SHA changed")
    seen={}
    ambiguous=set()
    qc=collections.Counter()
    for sheet,rows in original_rows(body):
        if not rows:continue
        head={str(v or '').strip().upper():k for k,v in rows[0].items()}
        required={"UKBMS CODE","EAST","NORTH","HABITAT1","HABITAT2"}
        if not required<=set(head):raise ValueError("Original crosswalk schema incomplete")
        for r in rows[1:]:
            qc["source_crosswalk_rows"]+=1
            site=key_site(r.get(head["UKBMS CODE"]))
            if not site:qc["unparseable_ukbms_code"]+=1;continue
            try:e=float(r.get(head["EAST"]));n=float(r.get(head["NORTH"]))
            except (TypeError,ValueError):qc["missing_east_north"]+=1;continue
            if not(math.isfinite(e) and math.isfinite(n) and 0<=e<=700000 and 0<=n<=1300000):
                qc["invalid_grid_coordinates"]+=1;continue
            habitat=tuple(str(r.get(head[k]) or "").strip().lower()[:80] for k in ("HABITAT1","HABITAT2"))
            if site in seen and (abs(seen[site][0]-e)>500 or abs(seen[site][1]-n)>500):
                ambiguous.add(site)
            else:seen[site]=(e,n,habitat)
    for s in ambiguous:seen.pop(s,None)
    qc["ambiguous_transect_coordinates"]=len(ambiguous)
    qc["original_unique_geocoded_crosswalk_transects"]=len(seen)
    return seen,dict(qc)

def source_visits():
    counts=collections.defaultdict(lambda:[0,0,set(),set()])
    source=set()
    for row in excel_rows(get_source("visits"),FIELDS_VISITS):
        event=parse_event(row)
        if event is None or event in source:continue
        source.add(event)
        site,year,month,day,t=event
        if year<1976 or year>2022 or not 4<=month<=9:continue
        x=counts[site]
        if year<=1995:
            x[0]+=1;x[2].add(year)
        if year>=2005:
            x[1]+=1;x[3].add(year)
    return counts

def evaluate(src,cross):
    strict={s for s,v in src.items() if v[0]>=10 and v[1]>=10 and len(v[2])>=3 and len(v[3])>=3}
    recurrent={s for s,v in src.items() if v[0]>0 and v[1]>0}
    match=set(cross)&set(src)
    geo_strict=strict & match
    habitat=collections.Counter()
    for s in geo_strict:
        for h in cross[s][2]:
            if h:habitat[h]+=1
    out={
      "schema":"chocho_uk_ltmn_original_transsect_site_geography_join_20261028_v01",
      "status":"SOURCE_CROSSWALK_AND_ORIGINAL_VISITS_PAIRED",
      "source_sha256":{"site_crosswalk":SHA},
      "original_unique_visitation_transects_1976_2022":len(src),
      "repeated_any_1976_95_vs_2005_22":len(recurrent),
      "strict_year_and_visit_matched_transects":len(strict),
      "source_geocoded_transects_matching_any_visits":len(match),
      "source_geocoded_transects_matching_strict_panel":len(geo_strict),
      "source_visits_on_strict_transects_1976_95":sum(src[s][0] for s in strict),
      "source_visits_on_strict_transects_2005_22":sum(src[s][1] for s in strict),
      "habitats_on_geocoded_strict_panel":dict(sorted(habitat.items())),
      "crosswalk_is_not_documented_contemporary_host_species":True,
      "zero_butterfly_count_is_not_inferred_from_missing_species_row":True,
      "original_site_numbers_easting_northing_and_habitat_by_site_exported":False,
    }
    return out

def run(outpath):
    c,q=read_crosswalk()
    src=source_visits()
    d=evaluate(src,c)
    d["source_crosswalk_qc"]=q
    outpath.parent.mkdir(parents=True,exist_ok=True)
    outpath.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(d,sort_keys=True))
    return d

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    run(p.parse_args().output)

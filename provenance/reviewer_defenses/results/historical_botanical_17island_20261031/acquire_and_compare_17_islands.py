#!/usr/bin/env python3
"""Extend source-pinned GBIF historical plants to ALL 17 nodes in 14 near island pairs.

Uses the original 23-port botanical 2001-2010 GBIF facet artifact (SHA locked),
fetches ONLY the 15 previously unobserved WPI-port/year histories and pools
indexed *speciesKey sets* across ports of each island. Indexed plants are not
host identity, biological absence, alien status, shipping, or establishment.
GBIF may be reindexed; using 2019 WPI locations is retrospectively geographic.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, math, sys, time
from collections import defaultdict, Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
B=Path("provenance/reviewer_defenses/results")
EARLY=B/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
PORT_STATUS=B/"port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
BASELINE=B/"maritime_host_steppingstone_20261030/near_island_record_order_baseline_20261030.json"
BUTTERFLIES=B/"island_host_network_20261030/island_host_network_pairs.csv"
ORIGINAL_AUDIT=B/"historical_host_controls_20261031/daily_botanical_all23ports_20261031.json"
ORIGINAL_SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
ORIGINAL_GIT_BLOBS={
 EARLY:"0a01831056c2cd59a4bd41bd1dbd4c84d77cf922",
 PORT_STATUS:"fdf73b11ddcb1644f46730123350ae52cc61cae6",
 BASELINE:"345c20c2c6492a29f5d23f06b5f95655f3e72f7c",
 BUTTERFLIES:"5663d2ffe42419f5c4496b92c9f3bbf76175907b",
 ORIGINAL_AUDIT:"4df7b05399a5c1324caffdd19c11f1553b592b9a",
}
ACQUIRE_MOD=B/"historical_host_controls_20261031/acquire_dated_botanical_facets.py"
NEW_COUNTRY_CODES={"CA":{"CA","CANADA"},"ES":{"ES","SPAIN"},"GR":{"GR","GREECE"}}
DATE_WINDOW=[2001,2010]
BASIS="2019 World Port Index port points, 5km GBIF geoDistance circular searches"
MIN_TAXA_FOR_BOTANICAL_COMPARISON=20
MIN_EARLY_YEARS_FOR_BOTANICAL_COMPARISON=2

def pinned(root,relative):
    data=(root/relative).read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if sha!=ORIGINAL_GIT_BLOBS[relative]:raise ValueError("Changed public source Git blob "+str(relative))
    return data.decode("utf-8")

def readcsv(s):
    return list(csv.DictReader(io.StringIO(s)))

def rank(v):
    idx=sorted(range(len(v)),key=lambda i:v[i]);z=[0.0]*len(v);k=0
    while k<len(idx):
        t=k+1
        while t<len(idx) and v[idx[t]]==v[idx[k]]:t+=1
        pos=(k+t-1)/2+1
        for i in idx[k:t]:z[i]=pos
        k=t
    return z

def cor(x,y):
    if len(x)!=len(y) or len(x)<3:return None
    mx=sum(x)/len(x);my=sum(y)/len(y)
    xx=sum((v-mx)**2 for v in x);yy=sum((v-my)**2 for v in y)
    if xx<=0 or yy<=0:return None
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(xx*yy)

def spearman(x,y):return cor(rank(x),rank(y))

def read_original(source_gzip):
    data=source_gzip.read_bytes()
    if hashlib.sha256(data).hexdigest()!=ORIGINAL_SHA:
        raise ValueError("Old 23-port plant source bytes changed")
    rows=readcsv(gzip.decompress(data).decode("utf-8"))
    if len(rows)!=9788:raise ValueError("Old 23-port plant species-year count changed")
    return [{ "island_id":str(x["island_id"]), "port_id":str(x["port_id"]),
        "year":int(x["year"]),"speciesKey":str(x["gbif_speciesKey"]),
        "indexed_records":int(x["indexed_record_count"])} for x in rows]

def acquire_new(root,remaining,out_gzip):
    import importlib.util
    spec=importlib.util.spec_from_file_location("source_2019_gbif_wpi",root/ACQUIRE_MOD)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.COUNTRIES.update(NEW_COUNTRY_CODES)
    codes={x["port_id"]:x["country"] for x in readcsv(pinned(root,PORT_STATUS))}
    ports=mod.get_ports(remaining,codes)
    if len(ports)!=15:raise ValueError("Did not resolve all 15 source WPI port identities")
    new=[];port_receipts=[]
    for p in sorted(remaining,key=lambda r:int(r["port_id"])):
        pid=int(p["port_id"]);iid=str(p["landmass_id"])
        source=ports[pid]
        counts=0;uniq=set();active=set()
        for y in range(2001,2011):
            n,keys,_=mod.year_facets(source,y)
            counts+=n;uniq.update(keys)
            if n:active.add(y)
            for k,v in keys.items():
                if int(v)<1:raise ValueError("GBIF source speciesKey facet invalid")
                new.append({"island_id":iid,"port_id":str(pid),"year":y,
                            "speciesKey":str(k),"indexed_records":int(v)})
            time.sleep(0.15)
        port_receipts.append({"island_id":iid,"port_id":pid,"country":source["expected_country"],
            "wpi_port_name":source["name"],
            "original_2001_2010_indexed_plant_records":int(p["early_plant_records"]),
            "new_annual_requery_indexed_records":counts,
            "new_indexed_species_keys":len(uniq),
            "observed_years":len(active),
            "old_and_requery_count_agree":counts==int(p["early_plant_records"])})
        print("NEW_15_PORT_YEAR",iid,pid,len(uniq),counts,flush=True)
    if len(port_receipts)!=15:raise ValueError("Historical 15-port acquisition incomplete")
    out_gzip.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(out_gzip,"wt",newline="",encoding="utf-8") as stream:
        w=csv.writer(stream)
        w.writerow(["island_id","port_id","year","speciesKey","indexed_records"])
        for x in new:w.writerow([x[k] for k in ("island_id","port_id","year","speciesKey","indexed_records")])
    return new,port_receipts

def compute(root,old,novel,port_receipts,new_source_sha):
    screen=json.loads(pinned(root,BASELINE))
    fixed=readcsv(pinned(root,BUTTERFLIES))
    earlyports=readcsv(pinned(root,EARLY))
    oldaudit=json.loads(pinned(root,ORIGINAL_AUDIT))
    if oldaudit["specieskey_year_artifact"]["sha256"]!=ORIGINAL_SHA:
        raise ValueError("Previously verified all23 plant source changed")
    links=screen["pair_summaries"]
    if len(links)!=14:raise ValueError("Original14-pair butterfly chronology changed")
    ids={str(i) for r in links for i in (r["island_a"],r["island_b"])}
    previous={str(i["island_id"]) for i in oldaudit["by_island_species_key_unions"]}
    if (len(ids),len(previous),len(ids-previous))!=(17,9,8):
        raise ValueError("Original17 node and 9 previously sampled island scope changed")
    newly_island={r["island_id"] for r in port_receipts}
    if newly_island!=ids-previous:raise ValueError("New port island IDs differ from 8 missing source islands")
    species=defaultdict(set);years=defaultdict(set);earliest=defaultdict(dict)
    for r in [*old,*novel]:
        i=r["island_id"];k=r["speciesKey"];y=int(r["year"])
        if i not in ids or not 2001<=y<=2010:raise ValueError("Unknown island or year")
        species[i].add(k);years[i].add(y)
        earliest[i][k]=min(earliest[i].get(k,y),y)
    if len(species)>17:raise ValueError("Excess island IDs")
    pairlookup={tuple(sorted((r["island_a"],r["island_b"]),key=int)):r for r in fixed}
    screenlookup={tuple(sorted((str(r["island_a"]),str(r["island_b"])),key=int)):r for r in links}
    result=[]
    for key,s in sorted(screenlookup.items(),key=lambda v:(int(v[0][0]),int(v[0][1]))):
        if key not in pairlookup:raise ValueError("Original butterfly 1275-pair source missing near pair")
        p=pairlookup[key];a,b=key
        A=species[a];C=species[b]
        na=len(A);nb=len(C);intersection=len(A&C);union=len(A|C)
        early_jac=intersection/union if union else None
        if abs(float(p["port_distance_km"])-s["port_representative_distance_km"])>1e-9:
            raise ValueError("Port representative geography changed")
        eligible=(na>=MIN_TAXA_FOR_BOTANICAL_COMPARISON and nb>=MIN_TAXA_FOR_BOTANICAL_COMPARISON and
                  len(years[a])>=MIN_EARLY_YEARS_FOR_BOTANICAL_COMPARISON and
                  len(years[b])>=MIN_EARLY_YEARS_FOR_BOTANICAL_COMPARISON)
        result.append({
            "island_a":a,"island_b":b,
            "port_representative_distance_km":round(float(p["port_distance_km"]),6),
            "plant_taxon_union_a":na,"plant_taxon_union_b":nb,
            "plant_recorded_years_a":len(years[a]),"plant_recorded_years_b":len(years[b]),
            "plant_taxon_overlap":intersection,"plant_taxon_union_pair":union,
            "historical_plant_jaccard":None if early_jac is None else round(early_jac,8),
            "prior_butterfly_jaccard":round(float(p["earlier_butterfly_jaccard"]),8),
            "new_butterfly_jaccard":round(float(p["new_butterfly_jaccard"]),8),
            "shared_new_butterfly_species_pairs":int(p["shared_new_species"]),
            "upstream_butterfly_rerecord_before_downstream_first":int(s["upstream_re_record_before_downstream_first_year"]),
            "botanical_reporting_quality_gate":bool(eligible),
            "plant_2001_2010_status":"NO_SPECIES_KEY_INDEX" if not union else "INDEXED_SPECIES_KEYS",
            "host_identity_status":"UNCLASSIFIED_EXCEPT_PREEXISTING_FIXED_FOCAL_LINKS",
        })
    if len(result)!=14 or sum(x["shared_new_butterfly_species_pairs"] for x in result)!=41:
        raise ValueError("Butterfly published pair population failed")
    full=[r for r in result if r["historical_plant_jaccard"] is not None]
    qualified=[r for r in full if r["botanical_reporting_quality_gate"]]
    def summary(rows):
        x=[r["historical_plant_jaccard"] for r in rows]
        n=[r["new_butterfly_jaccard"] for r in rows]
        old=[r["prior_butterfly_jaccard"] for r in rows]
        return {"pairs":len(rows),
            "rho_historical_plants_vs_new_butterfly_composition":spearman(x,n),
            "rho_historical_plants_vs_earlier_butterfly_composition":spearman(x,old),
            "source_not_independent_pairwise_samples":True}
    by_island=[{"island_id":i,"n_indexed_species_keys":len(species[i]),
        "years_with_indexed_botanical_records":len(years[i]),
        "n_original_source_ports":sum(p["landmass_id"]==i for p in earlyports),
        "source_provenance":"REUSED_PREVIOUS_23_PORTS" if i in previous else "NEW_15_PORT_YEAR_FACETS"}
        for i in sorted(ids,key=int)]
    # Similarity of pre-recorded plants, NOT butterfly host network similarity.
    output={"schema":"chocho_17island_historical_plant_record_composition_v1",
      "status":"EXPLORATORY_EARLIER_INDEXED_BOTANICAL_COMPOSITION_NOT_HOST_CAUSATION",
      "audit_date_utc":datetime.now(timezone.utc).isoformat(),
      "source_github_git_blobs":{str(k):v for k,v in ORIGINAL_GIT_BLOBS.items()},
      "private_23port_gzip_sha256":ORIGINAL_SHA,
      "new_15port_gzip_sha256":new_source_sha,
      "new_port_source_requeries":port_receipts,
      "source_plant_species_taxon_year_rows_new":len(novel),
      "source_plant_species_taxon_year_rows_original":len(old),
      "source_port_count_original":23,"source_port_count_new":15,
      "islands":by_island,"source_pairs":result,
      "descriptive_all_available":summary(full),"botanical_quality_gated":summary(qualified),
      "qualification_rule":{"min_species_keys_per_island":MIN_TAXA_FOR_BOTANICAL_COMPARISON,
          "min_years_with_species_key_per_island":MIN_EARLY_YEARS_FOR_BOTANICAL_COMPARISON,
          "fixed_before_inspecting_17island_outcomes":True},
      "no_exact_plant_alias_or_historical_native_introduced_status":True,
      "not_specific_larval_hosts":True,
      "not_host_absence_or_butterfly_establishment":True,
      "no_historical_port_vessel_routes_or_actual_displacement":True,
      "not_causal_relay_or_resource_rescue":True,
      "2019_wpi_port_proxy_and_changing_GBIF_index_warning":True
    }
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--original-23port-gzip",type=Path,required=True)
    p.add_argument("--new-15port-gzip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from pprint import pprint
    root=args.root
    base=json.loads(pinned(root,BASELINE))
    old=read_original(args.original_23port_gzip)
    early=readcsv(pinned(root,EARLY))
    wanted={str(p) for s in base["pair_summaries"] for p in (s["island_a"],s["island_b"])}
    audited=json.loads(pinned(root,ORIGINAL_AUDIT))
    oldids={str(p["island_id"]) for p in audited["by_island_species_key_unions"]}
    remaining=[p for p in early if p["landmass_id"] in wanted-oldids]
    if len(remaining)!=15:
        raise ValueError("Expected exactly 15 as-yet-unacquired historical WPI ports")
    novel, receipts=acquire_new(root,remaining,args.new_15port_gzip)
    out=compute(root,old,novel,receipts,hashlib.sha256(args.new_15port_gzip.read_bytes()).hexdigest())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("HISTORICAL_17_ISLAND_BOTANICAL",json.dumps({
        "near_pairs":len(out["source_pairs"]), "islands":len(out["islands"]),
        "new_15port_species_year_facets":len(novel),
        "all":out["descriptive_all_available"],
        "reporting_quality":out["botanical_quality_gated"],
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()

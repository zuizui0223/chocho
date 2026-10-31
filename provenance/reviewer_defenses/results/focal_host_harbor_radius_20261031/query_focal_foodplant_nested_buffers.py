#!/usr/bin/env python3
"""Source-verifiable GBIF historical foodplant radius audit on fixed WPI ports.

Re-query the *same* GBIF index at 1km and 5km on the exact public 2019 WPI
coordinates, 2001-2010 Plantae taxon speciesKey and all data-quality filters.
The nested-buffer differences are indexed flora location, NOT historic ship
introduction, local larval diet, plant abundance or butterfly colonization.
No individual GBIF occurrence, record ID, or protected host-network edges.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from collections import defaultdict
from pathlib import Path
from urllib.error import URLError,HTTPError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
PINS={
  "hosts":("provenance/reviewer_defenses/results/dated_focal_host_taxonomy_20261031/focal_historical_host_taxon_year_audit.json",
       "3048c769d250f1b2f8df97c45aa20c8b91bc573b"),
  "ports":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/fixed_geometric_sites.csv",
       "3373121189ac80b237d99e2a1b0a172101cfacb8"),
  "site_index":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_plant_site_facets.csv",
       "eb09a3f3ec053c00b43e8bc0678d6643ce02697b")
}
GBIF="https://api.gbif.org/v1/occurrence/search"
AGENT="chocho-focal-foodplant-1km-vs-5km-source-index-audit/20261031"

def read(key):
    path,expected=PINS[key]
    data=(ROOT/path).read_bytes()
    digest=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if digest!=expected:
        raise ValueError("Git source changed: "+path+" "+digest)
    if path.endswith(".json"):return json.loads(data)
    return list(csv.DictReader(io.StringIO(data.decode("utf8"))))

def request(url,params):
    address=url+"?"+urlencode(params)
    last=None
    for trial in range(4):
        try:
            req=Request(address,headers={"User-Agent":AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=40) as response:
                if response.status!=200:raise ValueError("GBIF HTTP="+str(response.status))
                result=json.load(response)
            if not isinstance(result,dict) or result.get("error"):
                raise ValueError("GBIF source result invalid/error "+str(result.get("error"))[:180])
            return result
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as err:
            last=err
            if isinstance(err,HTTPError) and err.code not in (408,425,429,500,502,503,504):
                break
            time.sleep(min(5,0.5*2**trial))
    raise ValueError("GBIF source retry exhausted "+repr(last))

def lookup_one(key,port,radius):
    params={"taxonKey":str(key),
        "year":"2001,2010","hasCoordinate":"true",
        "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
        "geoDistance":f"{port['latitude']},{port['longitude']},{radius}m",
        "facet":"speciesKey","facetLimit":"1000","limit":"0"}
    result=request(GBIF,params)
    count=int(result["count"])
    facets={str(x["field"]).upper():x.get("counts",[]) for x in result.get("facets",[])}
    if count and "SPECIES_KEY" not in facets:
        raise ValueError("Species-specific query returned positive records but no SPECIES_KEY")
    candidates={str(x.get("name")):int(x.get("count",0)) for x in facets.get("SPECIES_KEY",[])}
    # Should only contain the focal accepted species, since we explicitly use
    # the accepted species taxonKey, but save any GBIF resolution surprises.
    if len(candidates)>=1000:
        raise ValueError("Unexpected 1000+ descendant speciesKeys from species-level query")
    if sum(candidates.values())>count:raise ValueError("Indexed species facets exceed total")
    return {"indexed_taxonKey_records":count,
            "focal_accepted_speciesKey_indexed_records":int(candidates.get(str(key),0)),
            "all_facet_species_keys":candidates,
            "facet_unassigned_occurrences":count-sum(candidates.values())}

def run():
    hosts=read("hosts")["island_taxon_observations"]
    ports=read("ports")
    site=read("site_index")
    portmap={x["port_id"]:x for x in ports if x["site_type"]=="port"}
    site_map={x["site_id"]:x for x in site}
    assert len(portmap)==38
    combos=[]
    for p in hosts:
        for port_id in p["port_ids"]:
            port=portmap[str(port_id)]
            if port["landmass_id"]!=str(p["island_id"]):
                raise ValueError("Source focal host/port geographic mismatch")
            combos.append((p,port))
    assert len(combos)==11 and len(hosts)==4
    rows=[]
    for focal,port in combos:
        key=str(focal["gbif_species_key"])
        one=lookup_one(key,port,1000)
        five=lookup_one(key,port,5000)
        if one["indexed_taxonKey_records"]>five["indexed_taxonKey_records"]:
            raise ValueError("Unexpected non-nested GBIF radius counts")
        if one["focal_accepted_speciesKey_indexed_records"]>five["focal_accepted_speciesKey_indexed_records"]:
            raise ValueError("Non-nested focal accepted plant species")
        indexed=key in set(filter(None,site_map["PORT_"+port["port_id"]]["species_keys"].split(";")))
        # Comparison with the previously archived broad-kindgom 1km facet.
        # A mismatch is source-index drift, not an ecological change.
        concordance=(one["focal_accepted_speciesKey_indexed_records"]>0)==indexed
        rows.append({
            "host":focal["host_name"],"key":key,"island":focal["island_id"],
            "port_id":port["port_id"],
            "years":"2001-2010",
            "previous_5km_plant_indexed":True,
            "prior_5km_plant_first_index_year":focal["first_indexed_year"],
            "original_5km_plant_record_count_aggregate_island":focal["sum_indexed_records_across_potentially_overlapping_port_buffers"],
            "gbif_1km_taxonKey_records":one["indexed_taxonKey_records"],
            "gbif_1km_exact_accepted_speciesKey_records":one["focal_accepted_speciesKey_indexed_records"],
            "gbif_5km_taxonKey_records":five["indexed_taxonKey_records"],
            "gbif_5km_exact_accepted_speciesKey_records":five["focal_accepted_speciesKey_indexed_records"],
            "gbif_indexed_ring_1to5km_record_count_difference":
                 five["focal_accepted_speciesKey_indexed_records"]-one["focal_accepted_speciesKey_indexed_records"],
            "gbif_1km_two_source_query_concordant":concordance,
            "historic_plant_absence_inferred":False,
            "ship_transport_inferred":False
        })
        print("FOCAL_RADIUS",focal["host_name"],focal["island_id"],port["port_id"],
              one["focal_accepted_speciesKey_indexed_records"],
              five["focal_accepted_speciesKey_indexed_records"],flush=True)
        time.sleep(.07)
    if len(rows)!=11:raise ValueError("Focal source counts changed")
    by=defaultdict(list)
    for r in rows:by[(r["host"],r["island"])].append(r)
    summary=[{
        "host":name,"island":island,"port_sites":len(rs),
        "one_km_port_sites_with_index":sum(x["gbif_1km_exact_accepted_speciesKey_records"]>0 for x in rs),
        "five_km_port_sites_with_index":sum(x["gbif_5km_exact_accepted_speciesKey_records"]>0 for x in rs),
        "one_to_five_ring_sites_with_index":sum(x["gbif_indexed_ring_1to5km_record_count_difference"]>0 for x in rs),
        "all_nested_requeries_concordant_to_archived_1km":all(x["gbif_1km_two_source_query_concordant"] for x in rs)
        } for (name,island),rs in sorted(by.items())]
    receipt={
       "schema":"chocho_focal_butterfly_foodplant_nested_gbif_1km_5km_radius_v1",
       "status":"GBIF_SOURCE_REQUERY_NOT_CAUSAL_HARBOR_OR_HOST_RESCUE",
       "source_git_blobs":{k:{"path":path,"sha":sha} for k,(path,sha) in PINS.items()},
       "GBIF_API":GBIF,
       "species_specific_filter":"taxonKey=accepted plant species, validated by SPECIES_KEY facet",
       "window":"2001-2010",
       "quality_filters":["hasCoordinate=true","hasGeospatialIssue=false","occurrenceStatus=PRESENT"],
       "original_5km_source_focal_combinations":4,
       "focal_host_x_port_samples":len(rows),
       "host_island_summary":summary,
       "source_comparability_limits":[
          "1km and 5km queries are from same current GBIF index but may differ from the earlier first-indexed plant facet source.",
          "Multiple harbor buffers on one island may count the same plant observations, so raw record counts are not independent.",
          "Geospatial buffers are centered on modern 2019 WPI port coordinates and may include surroundings/urban vegetation.",
          "A plant 1-5 km from port is not evidence of ship movement, historical introduction, or local butterfly larval use.",
          "No GBIF record is biological proof of first establishment or true presence/absence at unindexed points.",
          "No unconditional butterfly species × island × year risk set or host-mediated movement endpoint is tested."
       ]
    }
    return rows,receipt

def main():
    rows,report=run()
    with (HERE/"focal_nested_harbor_queries.csv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    (HERE/"focal_nested_harbor_query_results.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":report["status"],
      "sample_count":report["focal_host_x_port_samples"],
      "host_island_summary":report["host_island_summary"]},indent=2,ensure_ascii=False))

if __name__=="__main__":main()

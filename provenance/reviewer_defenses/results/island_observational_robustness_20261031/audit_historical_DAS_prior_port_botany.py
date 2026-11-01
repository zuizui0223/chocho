#!/usr/bin/env python3
"""Historic VOC voyage endpoint evidence versus *indexed* early port plant data.

All ORIGINAL 51 WPI-port islands remain in the analytic frame. The
Huygens Dutch Asiatic Shipping 1595–1795 colonial voyage endpoint evidence
and GBIF 2001–2010 indexed port botanical observation opportunity are
SEPARATE data. Current candidate alien larval host and 2011–2022 indexed
butterfly richness are *descriptive geographic metadata*, not a biological
effect of seventeenth-century shipping.

A zero in 2001–2010 GBIF plant observations is UNOBSERVED plant records,
never botanical absence. No historical ship cargo/plant immigration/actual
butterfly breeding is identifiable from this cross-sectional join.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
SOURCE=BASE/"das_1595_1795_historic_island_endpoint_counts_v1.csv"
DAS=BASE/"das_1595_1795_historic_island_voyage_endpoints_v1.json"
PORTS=BASE.parent/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
WPI=BASE/"wpi_51island_194ports_geocoded.csv"
OUT=BASE/"das_hist_voyages_early_port_plant_record_coverage.json"


def rows(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):
            h.update(x)
    return h.hexdigest()


def audit():
    das=json.loads(DAS.read_text(encoding="utf8"))
    if das["schema"]!="chocho_source_pinned_DAS_Dutch_Europe_Asia_island_voyage_endpoint_v1":
        raise RuntimeError("Source original historic 8178-voyage receipt drift")
    islands=rows(SOURCE)
    plants=rows(PORTS)
    wpi=rows(WPI)
    if len(islands)!=51 or len(plants)!=194 or len(wpi)!=194:
        raise RuntimeError("Original historic island/194 indexed harbor records changed")
    ip={p["port_id"]:p["landmass_id"] for p in wpi}
    seen=set()
    early=defaultdict(lambda:{"records":0,"port_with_records":0,"ports":0})
    alltimes=defaultdict(int)
    for p in plants:
        port=p["port_id"]
        if port in seen or port not in ip or ip[port]!=p["landmass_id"]:
            raise RuntimeError("Nonunique or WPI-inconsistent early botanical port")
        seen.add(port)
        island=p["landmass_id"]
        earlier=int(p["early_plant_records"])
        latest=int(p["all_plant_records"])
        if earlier<0 or latest<earlier:
            raise RuntimeError("Early botanical records exceed entire source record scope")
        early[island]["records"]+=earlier
        early[island]["port_with_records"]+=int(earlier>0)
        early[island]["ports"]+=1
        alltimes[island]+=latest
    if sum(x["records"] for x in early.values())!=155002:
        raise RuntimeError("Original 194-port historic 2001-2010 plant observation opportunity checksum changed")
    if len(early)!=51 or len(seen)!=194:
        raise RuntimeError("Not every original 51 island has source-tested WPI port inventory")
    positive=[]
    absent=[]
    for r in islands:
        island=r["landmass_id"]
        voyage=int(r["historical_VOC_named_island_endpoint_total"])
        if voyage>0:
            j={
                "landmass_id":island,
                "historic_original_DAS_VOC_transoceanic_voyage_endpoint_rows":voyage,
                "historic_original_DAS_VOC_island_place_names":r["historical_source_port_labels"],
                "original_WPI_ports":early[island]["ports"],
                "2001_2010_indexed_plant_records_across_port_search_circles":early[island]["records"],
                "2001_2010_ports_with_positive_botanical_indexed_records":early[island]["port_with_records"],
                "2001_2024_indexed_candidate_alien_larval_host_species":int(r["port_indexed_candidate_alien_host_species_2001_2024"]),
                "2011_2022_newly_indexed_butterfly_species":int(r["newly_indexed_butterfly_species_2011_2022"]),
                "historic_NL_area_proxy_years":float(r["area_weighted_years_NL"]) if r["area_weighted_years_NL"] else None,
            }
            positive.append(j)
            if j["2001_2010_indexed_plant_records_across_port_search_circles"]==0:
                absent.append(island)
    if len(positive)!=4:
        raise RuntimeError("Precommitted city-level historic source support changed")
    return {
        "schema":"chocho_Huygens_DAS_voyage_endpoints_2001_2010_port_botanical_observation_support_v1",
        "status":"HISTORIC_SHIPPING_ENDPOINTS_OBSERVED_BUT_EARLY_GBIF_PLANT_OBSERVATION_NOT_COVERED_AT_ALL_HARBORS",
        "sources":{
            "historic_DAS_voyage_table_SHA256":sha(SOURCE),
            "historic_DAS_original_receipt_SHA256":sha(DAS),
            "port_GBIF_early_botanical_observation_counts_SHA256":sha(PORTS),
            "WPI_original_194_port_index_SHA256":sha(WPI),
        },
        "cohort":{"islands":51,"WPI_indexed_ports":194,
                  "source_botanical_observations_2001_2010":sum(x["records"] for x in early.values()),
                  "source_city_port_historic_VOC_endpoint_islands":len(positive),
                  "source_dated_historical_metropole_island_endpoint_rows":sum(x["historic_original_DAS_VOC_transoceanic_voyage_endpoint_rows"] for x in positive),
                  "source_historic_endpoint_islands_without_any_2001_2010_port_botanical_index":len(absent),
                  "source_historic_endpoint_islands_without_any_2001_2010_port_botanical_index_ids":absent},
        "matched_historical_source_islands":sorted(positive,key=lambda x:int(x["landmass_id"])),
        "limits":[
            "Historical VOC route names are matched to physical islands, NOT the exact 2019 WPI harbor locations.",
            "No original botanical GBIF records between 2001 and 2010 in a WPI 5km harbor buffer cannot establish plant absence during 17th century or 2001–2010.",
            "WPI buffer record totals can double count overlapping searches and may emphasize heavily sampled urban harbors.",
            "Modern host-plant alien status is not historic and the source does not contain voyage cargo species.",
            "Source islands were chosen for host/butterfly indexing in the initial positive observational study; ship linked islands are not independent replicates.",
            "This audit is source-coverage triage, NOT a historic trade effect on butterfly similarity.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    result=audit()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

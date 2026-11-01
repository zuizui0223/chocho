#!/usr/bin/env python3
"""Historical VOC ↔ modern AIS port ↔ prior GBIF plant coverage on frozen islands.

A DATED historic VOC journey endpoint (1595–1795), modern 2020 AIS
endpoint, and GBIF botanical recording opportunity (2001–2010) represent
DIFFERENT types of evidence. Their overlap by 51 physical island IDs is an
observation-coverage gate, not historic 2020 route persistence or proof that
host plants or butterflies were transported.

Original datasets/source and scientific selection remain frozen, with no
raw voyage, GBIF position, or botanical record made public here.
"""
from __future__ import annotations

from collections import Counter,defaultdict
import argparse
import csv
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
ISLANDS=BASE.parent/"island_host_network_20261030/island_host_network_islands.csv"
VOYAGE=BASE/"das_1595_1795_historic_island_endpoint_counts_v1.csv"
HIST_AUDIT=BASE/"das_1595_1795_historic_island_voyage_endpoints_v1.json"
PORTS=BASE/"wpi_51island_194ports_geocoded.csv"
BOTANY=BASE.parent/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
OUTPUT=BASE/"historical_VOC_modern_AIS_early_botany_51island_coverage_v1.json"
OUT_CSV=BASE/"historical_VOC_modern_AIS_early_botany_51island_coverage_v1.csv"


def rows(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit():
    nodes=[r for r in rows(ISLANDS) if r["eligible"]=="True"]
    if len(nodes)!=51:
        raise RuntimeError("Frozen 51 island source eligibility drift")
    ships={r["landmass_id"]:r for r in rows(VOYAGE)}
    meta=json.loads(HIST_AUDIT.read_text(encoding="utf8"))
    if (len(ships)!=51 or meta["source_original_8178_journey_count"]!=8178 or
        meta["metropole_to_island_voyage_endpoint_rows_total"]+
        meta["island_to_metropole_voyage_endpoint_rows_total"]!=6135 or
        meta["original_source_distinct_WPI51_island_pair_endpoint_routes"]):
        raise RuntimeError("Historic Huygens 8178 source grain changed")
    port_index={z["port_id"]:z["landmass_id"] for z in rows(PORTS)}
    bcount=Counter()
    bnport=Counter()
    botanical=rows(BOTANY)
    if len(port_index)!=194 or len(botanical)!=194:
        raise RuntimeError("WPI194 source has changed")
    seen=set()
    for row in botanical:
        pid=row["port_id"]
        if pid in seen or pid not in port_index or port_index[pid]!=row["landmass_id"]:
            raise RuntimeError("Original botanical WPI port coverage is not source-unique")
        seen.add(pid)
        bcount[row["landmass_id"]]+=int(row["early_plant_records"])
        bnport[row["landmass_id"]]+=int(int(row["early_plant_records"])>0)
    if sum(bcount.values())!=155002:
        raise RuntimeError("Early 2001–2010 port botanical record checksum drift")
    by_radius={}
    for km in (2,5,10):
        data=rows(BASE/f"ais_2020_wpi51_island_{km}km.WPI_port_AIS_source_coverage.csv")
        source={r["wpi_port_id"]:r for r in data}
        if len(source)!=194:
            raise RuntimeError("Original AIS/WPI endpoint cohort size mismatch")
        positive_ports=Counter()
        source_rows=Counter()
        for pid,landmass in port_index.items():
            if pid not in source or source[pid]["landmass_id"]!=landmass:
                raise RuntimeError("Original WPI/AIS port location source mismatch")
            score=int(source[pid]["ais2020_source_site_count"])
            cnt=int(source[pid]["ais2020_source_endpoint_rows"])
            if score>0:
                positive_ports[landmass]+=1
                source_rows[landmass]+=cnt
            if cnt>0 and score<=0:
                raise RuntimeError("Observed source endpoint without a matched site")
        by_radius[km]={"positive_ports":positive_ports,"source_rows":source_rows}
    records=[]
    for orig in sorted(nodes,key=lambda x:int(x["landmass_id"])):
        iid=orig["landmass_id"]
        hist=ships[iid]
        val={
            "landmass_id":iid,
            "original_WPI_ports":int(orig["ports"]),
            "historic_VOC_metropole_island_voyage_endpoints_1595_1795":
                int(hist["historical_VOC_named_island_endpoint_total"]),
            "historic_VOC_named_source_harbours":hist["historical_source_port_labels"],
            "area_weighted_Lenzner_four_empire_source_available":
                hist["has_four_empire_area_proxy"]=="1",
            "area_weighted_Lenzner_NL_duration_years":(
                float(hist["area_weighted_years_NL"])
                if hist["has_four_empire_area_proxy"]=="1" and hist["area_weighted_years_NL"]!="" else None
            ),
            "early_GBIF_port_botanical_index_rows_2001_2010":bcount[iid],
            "early_GBIF_port_buffer_with_at_least_one_plant_record":bnport[iid],
            "modern_catalogued_introduced_larval_host_species_2001_2024":
                int(orig["alien_host_species"]),
            "new_first_GBIF_indexed_butterfly_species_2011_2022":
                int(orig["new_butterfly_species"]),
        }
        for km in (2,5,10):
            val[f"AIS_2020_source_WPI_matched_ports_within_{km}km"]=(
                by_radius[km]["positive_ports"][iid])
            val[f"AIS_2020_source_endpoint_rows_within_{km}km"]=(
                by_radius[km]["source_rows"][iid])
        records.append(val)
    hist=[z for z in records
          if z["historic_VOC_metropole_island_voyage_endpoints_1595_1795"]>0]
    if len(hist)!=4:
        raise RuntimeError("Historic four city-name indexed islands changed")
    dis={}
    for km in (2,5,10):
        modern=[z for z in records if z[f"AIS_2020_source_WPI_matched_ports_within_{km}km"]>0]
        overlap=[z for z in hist if z[f"AIS_2020_source_WPI_matched_ports_within_{km}km"]>0]
        dis[str(km)]={
            "ais_source_island_nodes":len(modern),
            "historic_four_source_island_nodes_with_ais_endpoint":len(overlap),
            "historic_four_source_island_ids_with_ais_endpoint":[x["landmass_id"] for x in overlap],
            "historic_source_islands_with_ais_endpoint_AND_early_GBIF_plant_observation":
                len([x for x in overlap if x["early_GBIF_port_botanical_index_rows_2001_2010"]>0]),
            "triple_source_island_ids":[
                x["landmass_id"] for x in overlap
                if x["early_GBIF_port_botanical_index_rows_2001_2010"]>0],
        }
    report={
        "schema":"chocho_51island_historical_VOC_2020_AIS_2001_2010_botany_observational_coverage_v1",
        "status":"SOURCE_COHORT_OVERLAP_DIAGNOSIS_NOT_HISTORIC_ROUTE_OR_PLANT_EFFECT",
        "sources":{
            "island_51_sha256":sha(ISLANDS),
            "Huygens_8178_voyage_endpoint_island_summary_sha256":sha(VOYAGE),
            "original_harbour_WPI194_sha256":sha(PORTS),
            "2001_2010_plant_WPI194_sha256":sha(BOTANY),
            "AIS_2020_source_WPI194_2km_sha256":
                sha(BASE/"ais_2020_wpi51_island_2km.WPI_port_AIS_source_coverage.csv"),
            "AIS_2020_source_WPI194_5km_sha256":
                sha(BASE/"ais_2020_wpi51_island_5km.WPI_port_AIS_source_coverage.csv"),
            "AIS_2020_source_WPI194_10km_sha256":
                sha(BASE/"ais_2020_wpi51_island_10km.WPI_port_AIS_source_coverage.csv"),
        },
        "original_island_nodes":len(records),
        "original_2019_WPI_ports":len(port_index),
        "source_historic_DAS_city_endpoint_islands":len(hist),
        "source_historic_DAS_city_endpoint_island_ids":[x["landmass_id"] for x in hist],
        "historical_Dutch_voyage_source_islands_with_NL_area_proxy_25plus":[
            x["landmass_id"] for x in hist
            if x["area_weighted_Lenzner_NL_duration_years"] is not None
            and x["area_weighted_Lenzner_NL_duration_years"]>=25
        ],
        "historical_Dutch_voyage_source_islands_with_NL_area_proxy_zero":[
            x["landmass_id"] for x in hist
            if x["area_weighted_Lenzner_NL_duration_years"]==0
        ],
        "dated_historical_Dutch_voyage_endpoints_on_NL_area_proxy_zero_islands":sum(
            x["historic_VOC_metropole_island_voyage_endpoints_1595_1795"] for x in hist
            if x["area_weighted_Lenzner_NL_duration_years"]==0
        ),
        "historic_VOC_ports_with_2001_2010_GBIF_index_zero":len([
            x for x in hist if x["early_GBIF_port_botanical_index_rows_2001_2010"]==0]),
        "historic_city_endpoint_and_early_plant_observation_island_ids":[
            x["landmass_id"] for x in hist
            if x["early_GBIF_port_botanical_index_rows_2001_2010"]>0],
        "AIS_match_radius_sensitivity":dis,
        "hist_source_island_details":hist,
        "important_limits":[
            "Both Huygens historical voyage endpoints and 2020 AIS source endpoints are dated observations, not comparable link types or proof of route persistence.",
            "Lenzner's area-weighted Dutch imperial duration is different from documented Dutch VOC port exchange. Source Sri Lanka has historic Dutch VOC voyages but its selected botanical area proxy NL duration is zero; do not call this voyage-free or assume the occupation measure captures all trade.",
            "The Huygens 8178 voyage set has ZERO documented direct between-WPI51 island routes; it is not a complete historic intra-Asia port network.",
            "No pre2011 GBIF plant records in one port 5km source buffer means no indexed observation in that search, not actual absence of the plant or its cargo.",
            "2019 WPI port points matched 2020 AIS and 1595–1795 *island* place labels at different geographic grains.",
            "Historical port origin/destination source heavily concentrated in Java 6092/6135 voyages. An all-island correlation on four source-positive nodes would be biased.",
            "No examined record documents actual VOC-transported compatible larval host seeds, plants, eggs, larvae or butterfly reproduction.",
            "Source data are selected 51 islands with HOSTS/GBIF catalogues and 194 harbor WPI points; cannot generalize to all islands.",
        ]
    }
    return records,report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--report",type=Path,default=OUTPUT)
    p.add_argument("--region-csv",type=Path,default=OUT_CSV)
    args=p.parse_args()
    records,report=audit()
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with args.region_csv.open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps({k:v for k,v in report.items() if k not in ("sources","important_limits")},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

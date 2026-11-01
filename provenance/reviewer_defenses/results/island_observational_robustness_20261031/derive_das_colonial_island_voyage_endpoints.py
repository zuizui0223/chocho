#!/usr/bin/env python3
"""Published 1595–1795 DAS VOC journeys as historical island endpoint context.

Source: Huygens Dutch-Asiatic Shipping, normalized in exact IBM repository SHA
33ae4ef8dfa6acd5f35e619e079ceec518306c9f.
Different source from ~57,669 IISG intra-Asia voyage rows and 645,581
2020 AIS voyage rows. Never relabel transoceanic endpoint pairs as
intra-Asian direct shipping, cargo of plant propagules or entomological
invasion pathways.

Only historical SOURCE PORT VOCABULARY is matched to predeclared 51-island
physical landmasses; link to WPI island ID does not certify historical
operation of 2019 World Port Index harbor or colonial taxon arrival.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter,defaultdict
from datetime import date
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
SOURCE=Path("_historical_DAS/data/voyages.json")
LOCK=BASE/"das_island_historic_endpoint_lock_20261101.json"
HISTORY=BASE/"island_port_colonial_10km.islands.csv"
COHORT=BASE.parent/"island_host_network_20261030/island_host_network_islands.csv"
OUT=BASE/"das_1595_1795_historic_island_voyage_endpoints_v1.json"
OUT_ISLANDS=BASE/"das_1595_1795_historic_island_endpoint_counts_v1.csv"
SOURCE_SHA256="c875a0a59568713347b1d50e4a2f62708006be524133d30ee5ff4fe23f05b77f"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def year(value):
    s=str(value or "").strip()
    if len(s)<4 or not s[:4].isdigit():
        return None
    i=int(s[:4])
    return i if 1595<=i<=1795 else None


def get_source(source=SOURCE):
    if sha256(source)!=SOURCE_SHA256:
        raise RuntimeError("Normalized historic DAS source SHA drift")
    d=json.loads(source.read_text(encoding="utf8"))
    if not isinstance(d,list) or len(d)!=8178:
        raise RuntimeError("Source historic 8178 journeys missing or changed")
    if not all(isinstance(z,dict) for z in d):
        raise RuntimeError("Original source voyage rows malformed")
    return d


def validate_lock(lock=LOCK):
    d=json.loads(lock.read_text(encoding="utf8"))
    if d.get("schema")!="chocho_huygens_DAS_island_endpoints_prelocked_mapping_v1":
        raise RuntimeError("Historic island name mapping lock not source fixed")
    precise={k.casefold():v["island_id"] for k,v in d["historical_to_WPI51"]["precise"].items()}
    broad={k.casefold():v["island_id"] for k,v in d["historical_to_WPI51"]["coarse_island_labels_sensitivity_only"].items()}
    if set(precise)&set(broad):
        raise RuntimeError("Overlapping historical port aliases")
    home={name.casefold() for name in d["metropole_NL_historic_places"]}
    if (set(precise)|set(broad))&home:
        raise RuntimeError("Home harbor accidentally mapped to tropical island")
    return precise,broad,home


def classify(name,precise,broad,home,include_broad=False):
    label=str(name or "").strip().casefold()
    if label in precise:
        return "ISLAND_PRECISE",precise[label]
    if include_broad and label in broad:
        return "ISLAND_BROAD",broad[label]
    if label in home:
        return "DUTCH_HISTORIC_METROPOLE",None
    return "UNMAPPED_HISTORIC_PORT_OR_PLACE",None


def voyage_category(a,b):
    if a[0].startswith("ISLAND_") and b[0].startswith("ISLAND_"):
        return ("BETWEEN_WPI51_ISLAND_ENDPOINTS" if a[1]!=b[1]
                else "WITHIN_SAME_WPI51_ISLAND")
    if ((a[0]=="DUTCH_HISTORIC_METROPOLE" and b[0].startswith("ISLAND_")) or
        (b[0]=="DUTCH_HISTORIC_METROPOLE" and a[0].startswith("ISLAND_"))):
        return "DUTCH_METROPOLE_TO_OR_FROM_WPI51_ISLAND"
    if a[0].startswith("ISLAND_") or b[0].startswith("ISLAND_"):
        return "ISLAND_ENDPOINT_TO_UNMAPPED_OTHER_SOURCE_PORT"
    return "NO_WPI51_ISLAND_ENDPOINT"


def date_bin(year):
    if year<=1649:return "1595_1649"
    if year<=1699:return "1650_1699"
    if year<=1749:return "1700_1749"
    return "1750_1795"


def analyse(include_broad=False,source=SOURCE,lock=LOCK):
    precise,broad,home=validate_lock(lock)
    original=get_source(source)
    eligible={r["landmass_id"]:r for r in rows(COHORT) if r["eligible"]=="True"}
    history={r["landmass_id"]:r for r in rows(HISTORY)}
    if len(eligible)!=51 or len(history)!=51:
        raise RuntimeError("Frozen 51 WPI islands/historical source changed")
    for island in set(precise.values())|set(broad.values()):
        if island not in eligible:
            raise RuntimeError("Predeclared old VOC harbor lands outside original 51 islands: "+island)
    records=Counter()
    dated=Counter()
    distinct={k:set() for k in ("all","matched","arrival_confirmed")}
    historical_island_destinations=Counter()
    historical_island_origins=Counter()
    historic_pair_aggregate=Counter()
    port_labels=defaultdict(Counter)
    period_counts=defaultdict(Counter)
    arrival_status=Counter()
    for v in original:
        records["all_source_voyage_rows"]+=1
        t=year(v.get("departure_date"))
        a=classify(v.get("departure_port"),precise,broad,home,include_broad)
        b=classify(v.get("destination_port"),precise,broad,home,include_broad)
        cat=voyage_category(a,b)
        records[cat]+=1
        if t is None:
            records["no_valid_1595_1795_departure_year"]+=1
            continue
        distinct["all"].add(str(v.get("voyage_id","")))
        dated[cat]+=1
        arrival_time=year(v.get("arrival_date"))
        # A dated departure with a named destination does not prove arrival.
        arrived=(arrival_time is not None and t<=arrival_time<=t+3)
        if cat=="DUTCH_METROPOLE_TO_OR_FROM_WPI51_ISLAND":
            island=b[1] if a[0]=="DUTCH_HISTORIC_METROPOLE" else a[1]
            if a[0]=="DUTCH_HISTORIC_METROPOLE":
                historical_island_destinations[island]+=1
                direction="EUROPE_TO_ISLAND"
                port_labels[island][v["destination_port"]]+=1
            else:
                historical_island_origins[island]+=1
                direction="ISLAND_TO_EUROPE"
                port_labels[island][v["departure_port"]]+=1
            period_counts[island][date_bin(t)]+=1
            arrival_status["metropole_island_named_endpoint_rows"]+=1
            if arrived:
                arrival_status["metropole_island_arrival_dated_rows"]+=1
                distinct["arrival_confirmed"].add(str(v.get("voyage_id","")))
            else:
                arrival_status["metropole_island_arrival_unverified"]+=1
            distinct["matched"].add(str(v.get("voyage_id","")))
        if cat=="BETWEEN_WPI51_ISLAND_ENDPOINTS":
            pair=tuple(sorted((a[1],b[1]),key=int))
            historic_pair_aggregate[pair]+=1
    output_islands=[]
    for island in sorted(eligible,key=int):
        x=eligible[island]
        h=history[island]
        count_into=historical_island_destinations[island]
        count_out=historical_island_origins[island]
        row={
            "landmass_id":island,
            "port_indexed_candidate_alien_host_species_2001_2024":x["alien_host_species"],
            "early_butterfly_species_2001_2010":x["earlier_butterfly_species"],
            "newly_indexed_butterfly_species_2011_2022":x["new_butterfly_species"],
            "has_four_empire_area_proxy":h["eligible_for_history_test"],
            "area_weighted_years_NL":h["years_NL"] if h["eligible_for_history_test"]=="1" else "",
            "historical_VOC_named_island_endpoint_to":count_into,
            "historical_VOC_named_island_endpoint_from":count_out,
            "historical_VOC_named_island_endpoint_total":count_into+count_out,
            "historical_source_port_labels":";".join(sorted(port_labels[island])),
            "historic_1595_1649_voyage_rows":period_counts[island]["1595_1649"],
            "historic_1650_1699_voyage_rows":period_counts[island]["1650_1699"],
            "historic_1700_1749_voyage_rows":period_counts[island]["1700_1749"],
            "historic_1750_1795_voyage_rows":period_counts[island]["1750_1795"],
        }
        output_islands.append(row)
    positive=[r for r in output_islands if r["historical_VOC_named_island_endpoint_total"]>0]
    result={
        "schema":"chocho_source_pinned_DAS_Dutch_Europe_Asia_island_voyage_endpoint_v1",
        "status":"ORIGINAL_HUYGENS_VOC_HISTORICAL_EUROPE_ISLAND_VOYAGE_ENDPOINTS_NOT_DIRECT_INTERISLAND_LINKS",
        "source_huygens":"https://www.huygens.knaw.nl/en/resources/dutch-asiatic-shipping-in-the-17th-and-18th-centuries-2/",
        "source_normalized_IBM_commit":"33ae4ef8dfa6acd5f35e619e079ceec518306c9f",
        "source_sha256":sha256(source),
        "mapping_lock_sha256":sha256(lock),
        "source_original_8178_journey_count":len(original),
        "include_broad_geographic_island_labels":include_broad,
        "historical_source_rows_by_endpoint_type":dict(records),
        "historical_source_departure_dated_rows_by_endpoint_type":dict(dated),
        "distinct_historical_voyage_ids_all_with_dated_departure":len(distinct["all"]),
        "source_matched_metropole_to_island_arrival_date_checks":dict(arrival_status),
        "original_WPI_islands":len(output_islands),
        "islands_with_historical_Dutch_source_metropole_endpoint":len(positive),
        "historic_island_ids_with_dated_Europe_voyage_endpoint":[r["landmass_id"] for r in positive],
        "metropole_to_island_voyage_endpoint_rows_total":sum(historical_island_destinations.values()),
        "island_to_metropole_voyage_endpoint_rows_total":sum(historical_island_origins.values()),
        "original_source_distinct_WPI51_island_pair_endpoint_routes":{
            f"{a}|{b}":n for (a,b),n in sorted(historic_pair_aggregate.items(),key=lambda z:(int(z[0][0]),int(z[0][1])))},
        "original_source_Metropole_island_endpoint_voyage_ids":len(distinct["matched"]),
        "caveats":[
            "Huygens Dutch-Asiatic Shipping documents Europe–Asia *outward/homeward* voyages with Cape stopovers, not a complete seventeenth-century Asian coastal network.",
            "A source voyage destination label without dated arrival is not known to have arrived; separately audited.",
            "Only original 51-island physical landmass matches are shown; unmentioned islands do not have a certified absent historic route.",
            "Public historical port names are different from frozen 2019 WPI harbor identities. No port-specific historic WPI original ID is inferred.",
            "VOC vessel/cargo categories were not inspected for transported host plant species; no plant introduction or butterfly colonization is supported.",
            "The counts of 17th/18th century historical journeys cannot be used as the exposure to 2001–2022 butterfly indexing without multiple time and geographic controls.",
            "Contemporary host and butterfly richness on these islands are listed for transparent coverage only, not tested for causal shipping associations.",
        ],
    }
    return output_islands,result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--include-broad-island-labels",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    ap.add_argument("--island-csv",type=Path,default=OUT_ISLANDS)
    args=ap.parse_args()
    entries,report=analyse(args.include_broad_island_labels,args.source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with args.island_csv.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(entries[0]))
        w.writeheader();w.writerows(entries)
    print(json.dumps({k:v for k,v in report.items() if k!="caveats"},indent=2))


if __name__=="__main__":
    main()

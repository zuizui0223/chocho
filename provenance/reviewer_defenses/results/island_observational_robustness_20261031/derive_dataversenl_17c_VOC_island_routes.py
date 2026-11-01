#!/usr/bin/env python3
"""Direct historical VOC voyage record endpoints -> original 51 island landmasses.

Data: original 57,669 voyage rows, Leenstra & Parthesius 2022
doi:10.34894/5HGTCD SHA1 eed2479019a8570e392f37f3d781d5b0c67448b1.
The original 18 source field names are Dataverse .tab headers with DBF
type markers; from and to fields are VERPLAATS,C,52 and APLAATS,C,51.

Use ONLY a response-blind precommitted source city-port-to-landmass map.
Current port-indexed HOSTS/GBIF is consulted AFTER historical row aggregation.
No absent route is certified absent. Voyages are leg records, not cargo.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re

BASE=Path(__file__).resolve().parent
SOURCE=Path("_private_history/voc_57669_all_voyages.tab")
LOCK=BASE/"dataversenl_17c_to_WPI51_island_alias_lock_20261101.json"
SOURCE_LOCK=BASE/"dataversenl_17c_original_57669_lock_20261101.json"
ORIGINAL=BASE.parent/"island_host_network_20261030/island_host_network_pairs.csv"
ISLANDS=BASE.parent/"island_host_network_20261030/island_host_network_islands.csv"
EARLY_PLANTS=BASE.parent/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
OUT=BASE/"voc_57669_17c_interisland_routes_v1"
EXPECTED_SHA1="eed2479019a8570e392f37f3d781d5b0c67448b1"
EXPECTED_ROWS=57669
PORT_FROM="VERPLAATS,C,52"
PORT_TO="APLAATS,C,51"
DATE_FROM="VERTREK,D"
DATE_TO="AANKOMST,D"
SOURCE_KEY="REISN,N,6,0"


def digest(path,hashkind="sha256"):
    h=hashlib.new(hashkind)
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def source_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def normalize_port(name):
    return " ".join(str(name or "").strip().casefold().split())


def source_year(v):
    """Interpret the original VOC Dataverse dates, normally MM/DD/YYYY.

    Earlier parser incorrectly expected the year at the beginning,
    removing historical inter-island routes despite real source records.
    Reject impossible calendar days rather than guessing.
    """
    from datetime import datetime
    text = str(v or "").strip()
    if not text:
        return None
    if text.isdigit() and len(text)==4:
        y=int(text)
        return y if 1595<=y<=1795 else None
    for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%Y/%m/%d", "%Y%m%d"):
        try:
            y=datetime.strptime(text,fmt).year
            return y if 1595<=y<=1795 else None
        except ValueError:
            continue
    return None


def year_bin(y):
    if y<1650:return "1595_1649"
    if y<1700:return "1650_1699"
    if y<1750:return "1700_1749"
    return "1750_1795"


def load_aliases(lock=LOCK,include_broad=False):
    x=json.loads(Path(lock).read_text(encoding="utf8"))
    if x.get("schema")!="chocho_verified_17c_VOC_57669_place_alias_WPI51_landmass_20261101_v1":
        raise RuntimeError("Source-response precommitted original harbor-name lock missing")
    strict={normalize_port(k):v for k,v in x["primary_specific_named_harbor_aliases"].items()}
    broad={normalize_port(k):v for k,v in x["sensitivities_broad_island_or_area_names_not_ports"].items()}
    if set(strict)&set(broad):
        raise RuntimeError("A historical harbor was ambiguously labelled an entire island")
    return dict(strict,**(broad if include_broad else {})),strict,broad


def analyze(source=SOURCE,lock=LOCK,include_broad=False):
    if digest(source,"sha1")!=EXPECTED_SHA1:
        raise RuntimeError("Historical original voyage source SHA1 mismatch")
    if not SOURCE_LOCK.exists():
        raise RuntimeError("Original 57,669 public-source provenance lock absent")
    aliases,strict,broad=load_aliases(lock,include_broad)
    island_rows=source_rows(ISLANDS)
    eligible={r["landmass_id"] for r in island_rows if r["eligible"]=="True"}
    if len(eligible)!=51 or not set(aliases.values()) <= eligible:
        raise RuntimeError("Historic alias points beyond original 51 physical island nodes")
    pairrows=source_rows(ORIGINAL)
    if len(pairrows)!=1275:
        raise RuntimeError("Original butterfly island-pair source changed")
    island_lookup={r["landmass_id"]:r for r in island_rows}
    early_port_rows=source_rows(EARLY_PLANTS)
    if len(early_port_rows)!=194:
        raise RuntimeError("Original 194 WPI early plant-port inventory changed")
    plant_early=Counter()
    plant_positive_ports=Counter()
    for p in early_port_rows:
        plant_early[p["landmass_id"]]+=int(p["early_plant_records"])
        plant_positive_ports[p["landmass_id"]]+=int(int(p["early_plant_records"])>0)

    counts=Counter()
    edges=Counter()
    edge_ids=defaultdict(set)
    edge_year=Counter()
    directed_edges=Counter()
    island_endpoint=Counter()
    port_spellings=Counter()
    all_ids=set()
    matched_geography_ids=set()
    route_kinds=Counter()
    first_time=None
    last_time=None
    with Path(source).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        cols=set(reader.fieldnames or ())
        if not {PORT_FROM,PORT_TO,DATE_FROM,DATE_TO,SOURCE_KEY} <= cols or len(cols)!=18:
            raise RuntimeError("Archived historic voyage source schema drift")
        for r in reader:
            counts["source_record_rows"]+=1
            idval=str(r.get(SOURCE_KEY) or "").strip()
            if idval:
                all_ids.add(idval)
            fromname=normalize_port(r.get(PORT_FROM))
            toname=normalize_port(r.get(PORT_TO))
            a,b=aliases.get(fromname),aliases.get(toname)
            if a: island_endpoint[a]+=1
            if b: island_endpoint[b]+=1
            year=source_year(r.get(DATE_FROM))
            year_to=source_year(r.get(DATE_TO))
            if year is None:
                counts["source_departure_year_unparsed_or_outside_lock"]+=1
            else:
                first_time=year if first_time is None else min(year,first_time)
                last_time=year if last_time is None else max(year,last_time)
            if year is not None and year_to is not None and year_to<year:
                counts["arrival_year_precedes_departure_year"]+=1
            if a is None and b is None:
                route_kinds["NEITHER_ENDPOINT_MAPPED_TO_WPI51"]+=1
                continue
            if a is None or b is None:
                route_kinds["ONE_ENDPOINT_MAPPED"]+=1
                continue
            if a==b:
                route_kinds["SAME_PHYSICAL_ISLAND"]+=1
                continue
            route_kinds["TWO_DIFFERENT_WPI51_ISLANDS"]+=1
            if year is None:
                counts["two_island_route_rows_without_interpretable_date"]+=1
                continue
            edge=tuple(sorted((a,b),key=int))
            edges[edge]+=1
            directed_edges[(a,b)]+=1
            edge_year[(edge,year_bin(year))]+=1
            matched_geography_ids.add((a,b))
            if idval: edge_ids[edge].add(idval)
            port_spellings[(fromname,toname)]+=1
    if counts["source_record_rows"]!=EXPECTED_ROWS:
        raise RuntimeError("Original 57669 records not fully processed")
    if not edges:
        raise RuntimeError("No true 17th-century city-to-city island legs match prelocked aliases")
    positive_key=set(edges)
    original_lookup={
        tuple(sorted((r["island_a"],r["island_b"]),key=int)):r
        for r in pairrows
    }
    if len(original_lookup)!=1275:
        raise RuntimeError("Repeated old island pair source")
    endpoint_nodes={i for pair in edges for i in pair}
    pos_edges=[]
    pairs=[]
    for pair in sorted(original_lookup,key=lambda x:(int(x[0]),int(x[1]))):
        original=original_lookup[pair]
        n=edges.get(pair,0)
        data={
            "island_a":pair[0],"island_b":pair[1],
            "historic_source_17c_island_to_island_leg_positive":int(n>0),
            "historic_source_17c_voyage_leg_rows":str(n) if n else "",
            "historic_source_no_link_status":"NOT_OBSERVED_UNKNOWN" if not n else "DIRECT_LINK_OBSERVED",
            "port_indexed_candidate_host_jaccard":original["host_jaccard"],
            "GBIF_new_butterfly_jaccard_2011_2022":original["new_butterfly_jaccard"],
            "GBIF_earlier_butterfly_jaccard_2001_2010":original["earlier_butterfly_jaccard"],
            "source_port_distance_km":original["port_distance_km"],
        }
        pairs.append(data)
        if n:
            detail=dict(data)
            detail["historic_source_unique_voyage_record_ids"]=len(edge_ids[pair])
            detail["historic_period_1595_1649"]=edge_year.get((pair,"1595_1649"),0)
            detail["historic_period_1650_1699"]=edge_year.get((pair,"1650_1699"),0)
            detail["historic_period_1700_1749"]=edge_year.get((pair,"1700_1749"),0)
            detail["historic_period_1750_1795"]=edge_year.get((pair,"1750_1795"),0)
            detail["directed_legs_island_a_to_b"]=directed_edges.get(pair,0)
            detail["directed_legs_island_b_to_a"]=directed_edges.get((pair[1],pair[0]),0)
            detail["early_2001_2010_GBIF_plant_search_records_island_a"]=plant_early[pair[0]]
            detail["early_2001_2010_GBIF_plant_search_records_island_b"]=plant_early[pair[1]]
            detail["early_2001_2010_GBIF_port_searches_with_positive_plants_island_a"]=plant_positive_ports[pair[0]]
            detail["early_2001_2010_GBIF_port_searches_with_positive_plants_island_b"]=plant_positive_ports[pair[1]]
            pos_edges.append(detail)
    # Annotation of historical alias support is about the source rows only.
    counts["source_record_unique_REISN_ids"]=len(all_ids)
    n_alias_source=len(aliases)
    islands_covered=[]
    for iid in sorted(eligible,key=int):
        c=island_lookup[iid]
        islands_covered.append({
          "landmass_id":iid,
          "island_has_positive_historic_interisland_link":int(iid in endpoint_nodes),
          "historical_source_known_island_endpoint_rows":island_endpoint[iid],
          "island_17c_positive_source_neighbours":sum(iid in pair for pair in positive_key),
          "new_butterfly_indexed_species_2011_2022":int(c["new_butterfly_species"]),
          "port_indexed_candidate_alien_host_species":int(c["alien_host_species"]),
          "early_2001_2010_port_plant_search_records_total":plant_early[iid],
          "early_2001_2010_ports_with_plant_record":plant_positive_ports[iid],
        })
    count_positives=sum(r["historic_source_17c_island_to_island_leg_positive"] for r in pairs)
    if count_positives!=len(edges):
        raise RuntimeError("Historical island edge universe not fully mapped")
    report={
      "schema":"chocho_genuine_dataversenl_57669_17c_VOC_direct_island_edges_v1",
      "status":"SOURCE_VERIFIED_HISTORICAL_CITY_ENDPOINT_LINKS_NOT_HOST_PLANT_TRANSPORT_OR_BUTTERFLY_COLONIZATION",
      "original_source":{
        "doi":"10.34894/5HGTCD",
        "underlying_IISG":"hdl:10622/YP2COP",
        "historical_source_exact_SHA1":digest(source,"sha1"),
        "raw_rows_expected":EXPECTED_ROWS,
        "raw_rows_parsed":counts["source_record_rows"],
        "source_file_size_bytes":Path(source).stat().st_size,
        "precommitted_historic_city_to_landmass_alias_sha256":digest(lock),
        "original_host_butterfly_51_island_pair_sha256":digest(ORIGINAL),
      },
      "include_broad_island_level_names":include_broad,
      "prelocked_alias_count":n_alias_source,
      "alias_count_primary_specific_port":len(strict),
      "alias_count_broad_sensitivity_only":len(broad),
      "source_years_earliest_latest":[first_time,last_time],
      "source_quality":dict(counts),
      "source_endpoint_categories":dict(route_kinds),
      "historical_17c_positive_source_island_pairs":len(edges),
      "historical_17c_positive_source_directed_island_pairs":len(directed_edges),
      "historic_direct_island_pair_voyage_leg_rows_total":sum(edges.values()),
      "island_nodes_with_positive_17c_source_links":len(endpoint_nodes),
      "original_51_island_pairs":len(pairs),
      "original_51_island_nodes":len(eligible),
      "positive_direct_source_links_with_recorded_host_sharing":
          sum(float(r["port_indexed_candidate_host_jaccard"])>0 for r in pos_edges),
      "positive_direct_source_links_with_new_butterfly_sharing":
          sum(float(r["GBIF_new_butterfly_jaccard_2011_2022"])>0 for r in pos_edges),
      "positive_direct_source_links_with_old_butterfly_sharing":
          sum(float(r["GBIF_earlier_butterfly_jaccard_2001_2010"])>0 for r in pos_edges),
      "historically_linked_islands_with_no_2001_2010_port_plant_index":
          [r["landmass_id"] for r in islands_covered
           if r["island_has_positive_historic_interisland_link"] and
           r["early_2001_2010_port_plant_search_records_total"]==0],
      "source_named_historic_port_leg_aliases_top25":[
         {"from_port":a,"to_port":b,"rows":n}
         for (a,b),n in port_spellings.most_common(25)
      ],
      "limits":[
        "Published historical voyage leg records are not independent ships or cargo, REISN IDs may recur.",
        "The original voyage table is limited by Dutch VOC archival coverage and 17th-century sources.",
        "Current WPI index locates landmasses, not the very same harbor in 1600.",
        "Broad place-label sensitivity includes 'Ceylon','Taiwan',etc; these are not verified specific port names.",
        "Source-unmatched pair is UNKNOWN direct route status, not no voyage or zero transport.",
        "Early 2001–2010 port plant GBIF query count can be zero even for historically heavily connected islands; absence of record is not absence of plants.",
        "Current 2001–2024 host candidate lists may postdate 2011–2022 butterfly first GBIF indexing, and cannot date colonial introduction.",
        "No cargo plant identity, live larva/egg arrival, butterfly transport or reproduction is inferred.",
        "The matched link subset is concentrated in a few islands; pair-level p-values are not appropriate.",
      ],
    }
    return pos_edges,pairs,islands_covered,report


def save(prefix,positive,pairs,per_island,receipt):
    prefix.parent.mkdir(parents=True,exist_ok=True)
    for suffix,rows in ((".positive_historic_links.csv",positive),
                        (".all_51_island_pairs.csv",pairs),
                        (".island_source_coverage.csv",per_island)):
        path=Path(str(prefix)+suffix)
        with path.open("w",encoding="utf8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    Path(str(prefix)+".audit.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--broad-island-labels",action="store_true")
    ap.add_argument("--output-prefix",type=Path,default=OUT)
    args=ap.parse_args()
    a,b,c,report=analyze(args.source,LOCK,args.broad_island_labels)
    save(args.output_prefix,a,b,c,report)
    print(json.dumps({
        "status":report["status"],
        "archival_rows":report["source_quality"]["source_record_rows"],
        "source_endpoint_types":report["source_endpoint_categories"],
        "positive_island_routes":report["historical_17c_positive_source_island_pairs"],
        "positive_route_rows":report["historic_direct_island_pair_voyage_leg_rows_total"],
        "mapped_island_nodes":report["island_nodes_with_positive_17c_source_links"],
        "host_sharing_positive_routes":report["positive_direct_source_links_with_recorded_host_sharing"],
        "butterfly_sharing_positive_routes":report["positive_direct_source_links_with_new_butterfly_sharing"]
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

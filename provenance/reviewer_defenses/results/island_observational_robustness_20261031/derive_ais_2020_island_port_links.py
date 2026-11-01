#!/usr/bin/env python3
"""Exact public 2020 AIS-derived port-to-port voyage rows -> WPI island links.

Sources: Duvallet et al. (2024), Zenodo DOI 10.5281/zenodo.10783698,
MD5 af27e1bef1f706e2e29ed53aee47ff10, ZIP member triplet_2020.csv.

The CSV contains AIS-derived dated *from* and *to* harbor coordinates, which
are real OD *records in this dataset* (not complete global port-coverage).
Match those port coordinates to source-locked 194 World Port Index 2019 port
points on 51 islands using unique-nearest WGS84 geodesic criteria. A missing
route is UNKNOWN and NEVER a certified zero. Years are 2020, historically
after the original 2001–2010 butterfly baseline and some 2011–2022 events.

No original ship MMSI/IMO, dated voyage row, port name, or source raw
trajectory is committed. Derived anonymous 51-island positive route edges
are source-backed *present-day* associations only.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np

BASE=Path(__file__).resolve().parent
SOURCE=Path("_private_source_stage/triplet_2020.zip")
WPI=BASE/"wpi_51island_194ports_geocoded.csv"
PAIRS=BASE.parent/"island_host_network_20261030/island_host_network_pairs.csv"
ISLANDS=BASE.parent/"island_host_network_20261030/island_host_network_islands.csv"
LOCK=BASE/"ais_2020_island_voyage_lock_20261101.json"
OUTPUT_PREFIX=BASE/"ais_2020_wpi51_island"
EXPECTED_MD5="af27e1bef1f706e2e29ed53aee47ff10"
EXPECTED_CSV_BYTES=113764323
SOURCE_COLUMNS=(
    "date-out","date-in","mmsi_imo_vessel_mame",
    "time_mooring","nb_mooring_terminal",
    "Port_ID_waiting_area_id_from","PortName_from","lon_from","lat_from",
    "Port_ID_waiting_area_id_to","PortName_to","lon_to","lat_to",
)
R_EARTH_KM=6371.0088


def read_csv(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def md5sum(path):
    h=hashlib.md5()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def validate_wpi():
    ports=read_csv(WPI)
    if len(ports)!=194 or len({r["port_id"] for r in ports})!=194:
        raise RuntimeError("51-island 194 WPI port source changed")
    islands={r["landmass_id"] for r in read_csv(ISLANDS)
             if r["eligible"]=="True"}
    if len(islands)!=51 or {r["landmass_id"] for r in ports}!=islands:
        raise RuntimeError("AIS WPI port matching no longer covers frozen 51 landmass IDs")
    coords=np.array([[float(r["latitude"]),float(r["longitude"])] for r in ports])
    if np.any(~np.isfinite(coords)) or np.any(abs(coords[:,0])>90) or np.any(abs(coords[:,1])>180):
        raise RuntimeError("Original World Port Index coordinate not valid")
    return ports,coords


def read_ais(source=SOURCE):
    if not source.exists() or md5sum(source)!=EXPECTED_MD5:
        raise RuntimeError("Exact 2020 AIS-derived source ZIP MD5 not verified")
    with zipfile.ZipFile(source) as z:
        members=[x for x in z.infolist() if x.filename=="triplet_2020.csv"]
        if len(members)!=1 or members[0].file_size!=EXPECTED_CSV_BYTES:
            raise RuntimeError("2020 AIS source member changed")
        with z.open(members[0]) as b:
            f=io.TextIOWrapper(b,encoding="utf-8-sig",newline="",errors="replace")
            reader=csv.DictReader(f,delimiter=";")
            if not set(SOURCE_COLUMNS).issubset(set(reader.fieldnames or [])):
                raise RuntimeError("Expected 2020 AIS from-to port schema changed")
            sites={}
            route_count=Counter()
            diagnostics=Counter()
            endpoints=Counter()
            geo_lat=[90.,-90.]
            geo_lon=[180.,-180.]
            for r in reader:
                diagnostics["raw_source_route_rows"]+=1
                fromy=str(r.get("date-out") or "")[:4]
                toy=str(r.get("date-in") or "")[:4]
                if fromy!="2020" or toy!="2020":
                    diagnostics["excluded_cross_year_or_date_missing"]+=1
                    continue
                site=[]
                for side in ("from","to"):
                    sid=str(r.get("Port_ID_waiting_area_id_"+side) or "").strip()
                    try:
                        lat=float(r["lat_"+side])
                        lon=float(r["lon_"+side])
                    except (ValueError,TypeError):
                        break
                    if not (math.isfinite(lat) and math.isfinite(lon)
                            and abs(lat)<=90 and abs(lon)<=180):
                        break
                    name=str(r.get("PortName_"+side) or "").strip()
                    # Some source identifiers may label waiting areas, hence
                    # keep rounded coordinate with original source ID rather
                    # than combining distinct geographic points by name alone.
                    key=(sid,round(lon,3),round(lat,3))
                    sites[key]=(lon,lat)
                    site.append(key)
                    geo_lat[0]=min(geo_lat[0],lat)
                    geo_lat[1]=max(geo_lat[1],lat)
                    geo_lon[0]=min(geo_lon[0],lon)
                    geo_lon[1]=max(geo_lon[1],lon)
                if len(site)!=2:
                    diagnostics["excluded_missing_coordinate_or_invalid"]+=1
                    continue
                if site[0]==site[1]:
                    diagnostics["excluded_source_within_port_roundtrip"]+=1
                    continue
                route_count[(site[0],site[1])]+=1
                diagnostics["valid_2020_source_directed_voyage_rows"]+=1
                endpoints[site[0]]+=1
                endpoints[site[1]]+=1
            if diagnostics["raw_source_route_rows"]<1000:
                raise RuntimeError("Unexpectedly tiny AIS 2020 original source")
            return sites,route_count,dict(diagnostics),dict(endpoints),{
                "lat_min":geo_lat[0],"lat_max":geo_lat[1],
                "lon_min":geo_lon[0],"lon_max":geo_lon[1],
            }


def haversine_vector(lat,lon,orig):
    """Vectorized approximate WGS84 great circle distance in kilometers."""
    lat1=np.radians(float(lat))
    lon1=np.radians(float(lon))
    lat2=np.radians(orig[:,0])
    lon2=np.radians(orig[:,1])
    dl=lat2-lat1
    dn=lon2-lon1
    h=(np.sin(dl/2.)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dn/2.)**2)
    return 2*R_EARTH_KM*np.arcsin(np.sqrt(np.clip(h,0,1)))


def match_source_sites(sites,orig,threshold_km):
    if threshold_km not in (2,5,10):
        raise ValueError("Only predeclared 2,5,10km matches")
    result={}
    stats=Counter()
    selected_distance=[]
    ambiguous_nearest=[]
    for sid,(lon,lat) in sites.items():
        d=haversine_vector(lat,lon,orig)
        cand=np.flatnonzero(d<=threshold_km)
        if len(cand)==1:
            result[sid]=int(cand[0])
            selected_distance.append(float(d[int(cand[0])]))
            stats["unique_WPI_coordinate_match"]+=1
        elif len(cand)>1:
            stats["ambiguous_multiple_WPI_ports_in_radius"]+=1
            ambiguous_nearest.append(float(np.min(d)))
        else:
            stats["not_locatable_to_frozen_194_WPI_ports"]+=1
    return result,dict(stats),{
        "matched_source_location_max_km":max(selected_distance,default=None),
        "matched_source_location_median_km":float(np.median(selected_distance)) if selected_distance else None,
        "ambiguous_site_nearest_km_median":float(np.median(ambiguous_nearest)) if ambiguous_nearest else None,
    }


def to_derived_island_network(sites,routes,source_rows,port_matches,ports,pairs,threshold,source_zip_size):
    source_edges=Counter()
    island_edges=Counter()
    WPI_matches=Counter()
    source_matched_voyages=0
    source_same_island_voyages=0
    source_same_port_voyages=0
    for (s,t),count in routes.items():
        if s not in port_matches or t not in port_matches:
            continue
        wp1=ports[port_matches[s]]
        wp2=ports[port_matches[t]]
        WPI_matches[wp1["port_id"]]+=count
        WPI_matches[wp2["port_id"]]+=count
        source_matched_voyages+=count
        if wp1["port_id"]==wp2["port_id"]:
            source_same_port_voyages+=count
            continue
        if wp1["landmass_id"]==wp2["landmass_id"]:
            source_same_island_voyages+=count
            continue
        a,b=sorted((wp1["port_id"],wp2["port_id"]),key=int)
        source_edges[(a,b)]+=count
        a,b=sorted((wp1["landmass_id"],wp2["landmass_id"]),key=int)
        island_edges[(a,b)]+=count
    rows=[]
    for r in pairs:
        key=tuple(sorted((r["island_a"],r["island_b"]),key=int))
        c=island_edges.get(key,0)
        rows.append({
            "island_a":r["island_a"],"island_b":r["island_b"],
            "ais_2020_direct_voyages_observed_in_sample":c if c>0 else "",
            "has_positive_2020_AIS_observed_from_to_pair":int(c>0),
            "same_historical_pair_source":1,
            "island_2020_AIS_route_without_link_is_UNKNOWN":int(c==0),
        })
    if len(rows)!=1275:
        raise RuntimeError("Original 51-island pair cohort changed")
    a={
        "schema":"chocho_2020_AIS_WPI51_source_covered_positive_OD_v1",
        "status":"PUBLIC_AIS_OBSERVED_POSITIVE_VOYAGE_ROWS_NOT_COMPLETE_ISLAND_SHIPPING_CENSUS",
        "source":{
            "DOI":"10.5281/zenodo.10783698",
            "ZIP_MD5":EXPECTED_MD5,
            "archive_bytes":source_zip_size,
            "original_triplet_csv_bytes":EXPECTED_CSV_BYTES,
            "matched_original_WPI_port_locations":len(WPI_matches),
            "source_sites_total":len(sites),
        },
        "max_AIS_endpoint_to_single_WPI_km":threshold,
        "source_route_rows_audited":source_rows,
        "routes_with_both_WPI_endpoints_unique":source_matched_voyages,
        "matched_voyages_between_different_WPI_ports_on_same_island":source_same_island_voyages,
        "matched_voyages_between_different_islands":sum(island_edges.values()),
        "matched_directed_source_routes_distinct":len(routes),
        "source_undirected_WPI_port_pairs_observed":len(source_edges),
        "source_undirected_island_pairs_observed":len(island_edges),
        "original_island_pairs":1275,
        "original_island_nodes":51,
        "original_WPI_port_candidates":194,
        "island_nodes_with_observed_interisland_direct_voyages":len({
            x for ab in island_edges for x in ab}),
        "island_nodes_with_any_AIS_original_port_matches":len({
            ports[i]["landmass_id"] for i in set(port_matches.values())}),
        "island_pairs_with_no_observed_2020_AIS_voyage_UNCOVERED_OR_ZERO_UNKNOWN":
            len(rows)-len(island_edges),
        "source_data_limitations":[
            "AIS vessel message coverage differs by coast, ship type, provider, and time; missing edges are UNKNOWN, not no shipping.",
            "2020 direct journey rows were reconstructed by AIS authors; no independent trajectory-level revalidation or duplicate-vessel correction was done.",
            "From-to coordinates are matched uniquely within exact 2/5/10km radius to original 2019 WPI harbor points.",
            "A historical shipping route or a specific colonial-era botanical transfer was not measured.",
            "2020 voyages postdate some original 2011–2022 butterfly first-index events and 2001–2024 botanical indexing.",
            "All WPI endpoints of a matched AIS source row are mapped to original landmass IDs; no GBIF occurrence data are transferred.",
            "No observed route does not prove the route physically absent, including between two AIS-covered ports.",
        ],
    }
    portedges=[{
        "wpi_port_id_a":a,"wpi_port_id_b":b,
        "ais_2020_sampled_voyage_rows":c,
    } for (a,b),c in sorted(source_edges.items(),key=lambda kv:(int(kv[0][0]),int(kv[0][1])))]
    return rows,portedges,a


def write_csv(path,records):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf8",newline="") as f:
        if not records:
            raise ValueError("No qualifying positive route rows to export")
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--zip",type=Path,default=SOURCE)
    ap.add_argument("--output-prefix",type=Path,default=OUTPUT_PREFIX)
    args=ap.parse_args()
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock.get("schema")!="chocho_2020_AIS_source_od_WPI194_island51_v1":
        raise RuntimeError("Pre-analysis AIS vessel port matching source lock absent")
    sites,routecounts,rawcounts,endpoints,bbox=read_ais(args.zip)
    ports,coords=validate_wpi()
    sourcepairs=read_csv(PAIRS)
    summaries={}
    for distance in (2,5,10):
        mapping,matchstats,matchdist=match_source_sites(sites,coords,distance)
        rows,portedges,a=to_derived_island_network(
            sites,routecounts,rawcounts,mapping,ports,sourcepairs,distance,
            args.zip.stat().st_size)
        a["source_WPI_matches"]=matchstats
        a["source_geographic_bbox"]=bbox
        a["site_match_distance_diagnostics"]=matchdist
        a["AIS_original_source_unique_port_site_locations"]=len(sites)
        a["AIS_validated_original_voyage_edges"]=len(routecounts)
        # Freeze which original WPI port locations were actually present in
        # this AIS source. A non-match is UNOBSERVED, NOT "no harbor route".
        wpi_sites=defaultdict(list)
        for ais_site,wpi_idx in mapping.items():
            wpi_sites[wpi_idx].append(ais_site)
        coverage=[]
        for idx,p in enumerate(ports):
            matched_sites=wpi_sites.get(idx,[])
            coverage.append({
                "wpi_port_id":p["port_id"],
                "landmass_id":p["landmass_id"],
                "country":p["country"],
                "ais2020_source_site_count":len(matched_sites),
                "ais2020_source_endpoint_rows":sum(endpoints.get(x,0) for x in matched_sites),
                "ais2020_endpoint_source_status":(
                    "AIS_ENDPOINT_PRESENT_IN_SOURCE" if matched_sites
                    else "NOT_OBSERVED_OR_NOT_GEOCODE_MATCHED"),
            })
        a["island_nodes_with_positive_AIS_source_endpoint"]=len({
            v["landmass_id"] for v in coverage
            if v["ais2020_source_site_count"]>0})
        a["matched_source_port_sites_source_endpoint_rows"]=sum(
            v["ais2020_source_endpoint_rows"] for v in coverage)
        a["WPI_geocoded_source_sha256"]=sha256(WPI)
        a["island_pair_source_sha256"]=sha256(PAIRS)
        a["preanalysis_lock_sha256"]=sha256(LOCK)
        prefix=Path(str(args.output_prefix)+f"_{distance}km")
        write_csv(Path(str(prefix)+".island_pairs.csv"),rows)
        write_csv(Path(str(prefix)+".WPI_port_AIS_source_coverage.csv"),coverage)
        if portedges:
            write_csv(Path(str(prefix)+".observed_port_links.csv"),portedges)
        Path(str(prefix)+".audit.json").write_text(
            json.dumps(a,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
        summaries[str(distance)]={
            "matched_source_sites":len(mapping),
            "matched_WPI_ports":a["source"]["matched_original_WPI_port_locations"],
            "matched_island_nodes":a["island_nodes_with_any_AIS_original_port_matches"],
            "positive_island_pairs":a["source_undirected_island_pairs_observed"],
            "positive_interisland_voyage_rows":a["matched_voyages_between_different_islands"],
        }
    print(json.dumps({
        "source_year":"2020","source_raw_route_counts":rawcounts,
        "unique_source_port_sites":len(sites),
        "unique_directed_source_routes":len(routecounts),
        "matched_by_radius":summaries
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Map 2019 WPI island ports to Lenzner 2022 historical empire-duration proxy.

All sources pinned before outcome tests. Original 51 island × 194 WPI port
coordinates, historical area-weighted GloNAF->TDWG 369 WGSRPD3 durations
(GB, ESP, PRT, NL), exact official level3 boundaries.

2019 ports are NOT proof of an operating harbor in colonial times.
Area-weighted historic empire occupation of mapped botanical regions is NOT
a port-specific shipping record, plant introduction or butterfly arrival.
Never impute missing historical polygon coverage as no colonial occupation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from pyproj import Geod
from shapely.geometry import Point, shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree
from shapely.validation import make_valid

BASE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[4]
PORTS=BASE/"wpi_51island_194ports_geocoded.csv"
ISLANDS=BASE.parent/"island_host_network_20261030/island_host_network_islands.csv"
HISTORY=Path("_external/colonial/provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/lenzner_wgsrpd3_colonial_area_crosswalk.csv")
GEOMETRY=Path("_external/tdwg/geojson/level3.geojson")
OFFICIAL=ROOT/"data/frozen/wgsrpd3_level1_lookup_v0.1.json"
LOCK=BASE/"island_port_colonial_design_lock_20261101.json"
YEARS=("years_GB","years_ESP","years_PRT","years_NL")
GEOD=Geod(ellps="WGS84")


def rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def checksum(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def build_geometry(path=GEOMETRY, official_path=OFFICIAL):
    official=json.loads(official_path.read_text(encoding="utf8"))
    if official.get("source_commit")!="52da7828aba9d461dd133c27b3bd7a4407161f54":
        raise ValueError("TDWG official source pin mismatch")
    codes=set(official["units"])
    if len(codes)!=369:
        raise ValueError("Wrong frozen official regions count")
    source=json.loads(path.read_text(encoding="utf8"))
    members=[]
    repaired={}
    for record in source["features"]:
        name=str(record["properties"].get("LEVEL3_COD") or "").strip()
        if not name:
            continue
        poly=shape(record["geometry"])
        if poly.is_empty:
            raise RuntimeError("Official polygon empty for "+name)
        if not poly.is_valid:
            origarea=float(poly.area)
            poly=make_valid(poly)
            if poly.is_empty or not poly.is_valid or poly.area<=0:
                raise RuntimeError("Cannot source-repair official polygon "+name)
            repaired[name]=abs(poly.area-origarea)/max(origarea,1e-12)
        members.append((name,poly))
    members.sort(key=lambda x:x[0])
    ordered=[x[0] for x in members]
    if len(members)!=369 or set(ordered)!=codes or len(set(ordered))!=369:
        raise ValueError("Original official polygons not exactly frozen 369 units")
    polygons=[x[1] for x in members]
    return STRtree(polygons),polygons,ordered,repaired


def closest_tdgw_polygon(lon,lat,tree,polygons,codes,max_distance=25.):
    """Direct point coverage or nearest coast polygon within a fixed tolerance.

    Candidate search first queries within two geographic degrees. True nearest
    coastal distances are WGS84 Geod; source-resolution edge artifacts can
    move a 2019 WPI point slightly out of land. No >25km forced matching.
    """
    if not (-180<=lon<=180 and -90<=lat<=90):
        raise ValueError("WPI port coordinate is invalid")
    pt=Point(lon,lat)
    idxs=tree.query(pt,predicate="intersects")
    if len(idxs):
        names=sorted(codes[int(i)] for i in idxs)
        return names[0],0.,len(names),"DIRECT_COVERS"
    candidates=[int(i) for i in tree.query(pt.buffer(2.0))]
    if not candidates:
        candidates=[int(tree.nearest(pt))]
    distances=[]
    for i in set(candidates):
        try:
            nearest=nearest_points(pt,polygons[i])[1]
            _,_,distance=GEOD.inv(lon,lat,nearest.x,nearest.y)
        except (ValueError,TypeError):
            continue
        if math.isfinite(distance):
            distances.append((distance/1000.,codes[i]))
    if not distances:
        return None,None,0,"NO_POLYGON_NEARBY"
    distances.sort()
    best,n=distances[0]
    if best>max_distance:
        return None,best,0,"PORT_OFFSHORE_BEYOND_25KM"
    return n,best,1,"COASTAL_NEAREST"


def build(max_port_distance_km=10., quality_lower=.8, quality_upper=1.05,
          ports=PORTS, islands=ISLANDS, history=HISTORY, geometry=GEOMETRY,
          lookup=OFFICIAL):
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock["schema"]!="chocho_51island_194WPIports_x_Lenzner4empires_x_2010bilateral_pretest_v1":
        raise ValueError("Pre-response history test lock absent")
    if max_port_distance_km not in (2.,10.,25.):
        raise ValueError("Matched geocoding thresholds must be predeclared")
    tree,polygons,codes,repaired=build_geometry(geometry,lookup)
    allports=rows(ports)
    allislands={r["landmass_id"]:r for r in rows(islands) if r["eligible"]=="True"}
    historic={r["wgsrpd3_code"]:r for r in rows(history)}
    if len(allports)!=194 or len(allislands)!=51 or len(historic)!=369:
        raise RuntimeError("Existing frozen 51-island, 194-port or 369-region source changed")
    if set(p["landmass_id"] for p in allports)!=set(allislands):
        raise RuntimeError("Original ports do not match all 51 eligible landmasses")
    port_counts=Counter(p["landmass_id"] for p in allports)
    for landmass, source in allislands.items():
        if int(source["ports"])!=port_counts[landmass]:
            raise RuntimeError("Number of indexed WPI ports changed for island "+landmass)
    portrows=[]
    per_island=defaultdict(list)
    diagnostics=Counter()
    for p in allports:
        lid=p["landmass_id"]
        region,distance,n_direct,status=closest_tdgw_polygon(
            float(p["longitude"]),float(p["latitude"]),tree,polygons,codes)
        accepted=bool(region and distance<=max_port_distance_km)
        good=False
        cov=None
        if accepted:
            h=historic[region]
            cov=float(h["glonaf_covered_area_fraction_unclipped"])
            good=(quality_lower<=cov<=quality_upper)
        if not accepted:
            state="NOT_MAPPED_WITHIN_CUTOFF"
        elif not good:
            state="GLO_NAF_AREA_COVERAGE_NOT_ACCEPTABLE"
        else:
            state="SOURCE_AREA_QUALITY_PASS"
        r={
            "port_id":p["port_id"],"landmass_id":lid,
            "source_current_country":p["country"],"wgsrpd3_code":region or "",
            "port_to_polygon_geodesic_km":distance if distance is not None else "",
            "port_geometry_method":status,"multiple_direct_regions":n_direct,
            "source_GLONAF_area_coverage":cov if cov is not None else "",
            "colonial_source_acceptance":state,
        }
        for y in YEARS:
            if good:
                val=float(historic[region][y])
                if not math.isfinite(val) or val<0:
                    raise RuntimeError("Bad historic empire duration for "+region+"/"+y)
                r[y]=val
            else:
                r[y]=""
        portrows.append(r)
        per_island[lid].append(r)
        diagnostics[state]+=1
    islandrows=[]
    for lid in sorted(allislands,key=int):
        p=per_island[lid]
        accepted=[r for r in p if r["colonial_source_acceptance"]=="SOURCE_AREA_QUALITY_PASS"]
        countries=sorted({q["source_current_country"] for q in p})
        used=len(accepted)
        fraction=used/len(p)
        good=used>0 and fraction>=.5
        row={
            "landmass_id":lid,"WPI_ports_original":len(p),
            "WPI_ports_colonial_mapped":used,
            "WPI_port_history_coverage_fraction":fraction,
            "eligible_for_history_test":int(good),
            "current_WPI_country_codes":";".join(countries),
            "multi_country_island":int(len(countries)>1),
            "region_codes_mapped":";".join(sorted({r["wgsrpd3_code"] for r in accepted})),
            "mapped_WPI_unique_region_codes":len({r["wgsrpd3_code"] for r in accepted}),
            "max_mapped_coastal_distance_km":(
                max(float(r["port_to_polygon_geodesic_km"]) for r in accepted)
                if used else ""),
        }
        for empire in YEARS:
            if good:
                values=[float(r[empire]) for r in accepted]
                row[empire]=sum(values)/len(values)
                row[empire+"_minimum"]=min(values)
                row[empire+"_maximum"]=max(values)
            else:
                row[empire]=""
                row[empire+"_minimum"]=""
                row[empire+"_maximum"]=""
        islandrows.append(row)
    eligible=[x for x in islandrows if x["eligible_for_history_test"]==1]
    summary={
        "schema":"chocho_WPI_194ports_to_colonial_51islands_audited_v1",
        "status":"GEOCODED_EMPIRE_REGION_HISTORY_NOT_HISTORIC_PORT_SHIPPING",
        "pretest_lock_sha256":checksum(LOCK),
        "source_sha256":{
            "WPI_194ports_geocoded":checksum(ports),
            "original_51_island_source":checksum(islands),
            "Lenzner_369_area_weighted_empire_years":checksum(history),
            "TDWG_369_original_geojson":checksum(geometry),
        },
        "source_WPI_ports":len(allports),"source_islands":len(allislands),
        "WPI_to_GloNAF_port_status":dict(diagnostics),
        "original_geographic_regions":len(codes),
        "repaired_source_polygon_count":len(repaired),
        "repaired_source_polygon_ids":sorted(repaired),
        "per_island_matching_threshold_km":max_port_distance_km,
        "required_source_GLONAF_area_coverage":[quality_lower,quality_upper],
        "eligible_historic_islands":len(eligible),
        "island_ids_with_eligible_source":sorted([v["landmass_id"] for v in eligible],key=int),
        "mapped_WPI_ports_in_eligible_islands":sum(r["WPI_ports_colonial_mapped"] for r in eligible),
        "source_ambiguous_multi_country_islands":[
            v["landmass_id"] for v in islandrows if v["multi_country_island"]],
        "limitations":[
            "WPI coordinates catalogued near 2019, not verified harbors during historic occupation.",
            "Lenzner source region colonial duration is an area-weighted botanical-geography proxy, not port-specific economic or shipping exchange.",
            "Source GloNAF regions and TDWG regions have area-overlap mismatch; low-coverage ports excluded.",
            "Port coordinates can fall near a neighboring peninsula/shore; nearest TDWG region does NOT prove actual landmass identity.",
            "Currently missing colonial source support is UNKNOWN, not zero historical occupation.",
            "Island histories are averaged over original matched ports with equal weight, not port freight or passenger-volume weight.",
            "50% source-port coverage and 10km distance fixed independently of butterfly/host pair outcomes.",
            "Countries coded current as a location attribute; not historic colonial empire assignment.",
        ],
    }
    return portrows,islandrows,summary


def write_csv(path,entries):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(entries[0]))
        w.writeheader()
        w.writerows(entries)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--max-port-to-region-km",type=float,default=10.)
    p.add_argument("--ports",type=Path,default=PORTS)
    p.add_argument("--islands",type=Path,default=ISLANDS)
    p.add_argument("--history",type=Path,default=HISTORY)
    p.add_argument("--geometry",type=Path,default=GEOMETRY)
    p.add_argument("--official",type=Path,default=OFFICIAL)
    p.add_argument("--output-prefix",type=Path,default=BASE/"island_port_colonial_10km")
    args=p.parse_args()
    pp,ii,a=build(args.max_port_to_region_km,ports=args.ports,
                  islands=args.islands,history=args.history,
                  geometry=args.geometry,lookup=args.official)
    write_csv(Path(str(args.output_prefix)+".ports.csv"),pp)
    write_csv(Path(str(args.output_prefix)+".islands.csv"),ii)
    Path(str(args.output_prefix)+".audit.json").write_text(
        json.dumps(a,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in a.items() if k!="source_sha256"},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

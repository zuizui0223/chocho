#!/usr/bin/env python3
"""Map all available WPI port points to source-locked global GSHHG coastlines.

This is a descriptive geographic index of *modern* WPI facilities, not ship
traffic intensity, historical port introduction, butterfly arrival, or larval
foodplant supply. Island sample is preselected from >=20 km² GSHHG L1 shapes,
independent of later butterfly index events.

Use exact near-coast geodesic distances to original GSHHG polygons to avoid
mistaken assignments to nearby mainland or to other islands. Never infer
"port absent" where WPI has no index; that is only WPI unindexed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Point,box
from shapely.ops import nearest_points
from shapely.strtree import STRtree
from pyproj import Geod

from build_gshhg_global_island_universe import source_shape,GSHHG_SHA,AREA_MIN,AREA_MAX

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
FRAME=HERE/"global_island_frame_area20.csv"
FRAME_BLOB=None
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
USER_AGENT="chocho-20261031-GSHHG-global-WPI-modern-ports/1.0"
GEOD=Geod(ellps="WGS84")

def get(url,params):
    target=url+"?"+urlencode(params,doseq=True)
    err=None
    for it in range(4):
        try:
            with urlopen(Request(target,headers={"User-Agent":USER_AGENT,
                     "Accept":"application/json"}),timeout=50) as f:
                result=json.load(f)
            if not isinstance(result,dict) or result.get("error"):
                raise ValueError("Invalid WPI source JSON "+repr(result.get("error"))[:300])
            return result
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as ex:
            err=ex
            if isinstance(ex,HTTPError) and ex.code not in (408,425,429,500,502,503,504):
                break
            time.sleep(min(5,.5*2**it))
    raise RuntimeError("WPI source failed "+str(err)[:350])

def all_wpi():
    offset=0;rows=[]
    while True:
        q=get(WPI,{"where":"1=1","outFields":"OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
                   "returnGeometry":"false","f":"json","resultOffset":offset,
                   "resultRecordCount":1000,"orderByFields":"OBJECTID ASC"})
        part=q.get("features",[])
        if not part:break
        rows.extend(x["attributes"] for x in part)
        offset+=len(part)
        if offset>20000:raise ValueError("Unexpected all WPI source size")
        if len(part)<1000 and not q.get("exceededTransferLimit"):break
    p=[]
    for r in rows:
        try:
            key=int(r["OBJECTID"]);lat=float(r["LATITUDE"]);lon=float(r["LONGITUDE"])
        except (TypeError,KeyError,ValueError):
            continue
        if not (-90<=lat<=90 and -180<=lon<=180):continue
        p.append({"port_id":key,"port_name":str(r.get("PORT_NAME") or ""),
                  "country":str(r.get("COUNTRY") or ""),"longitude":lon,"latitude":lat})
    if len(p)<3000 or len({x["port_id"] for x in p})!=len(p):
        raise ValueError("World Port Index catalog has invalid or duplicate original WPI ids")
    return sorted(p,key=lambda x:x["port_id"]),len(rows)

def read_frame():
    if not FRAME.is_file():
        raise ValueError("Prior globally outcome-independent island frame not generated")
    raw=FRAME.read_bytes()
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf8"))))
    ids={r["landmass_id"] for r in rows}
    if len(ids)<3000 or len(ids)!=len(rows):
        raise ValueError("Expected >3000 unique global GSHHG candidate islands")
    return rows,ids,hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def nearest_original_landmass(port,tree,geoms,land_ids):
    lon=float(port["longitude"]);lat=float(port["latitude"])
    point=Point(lon,lat)
    # Covers >5km in geographic coordinate box, including polar WPI ports.
    rlat=.06
    rlon=min(5.0,max(.06,5.8/(111.1*max(.008,abs(math.cos(math.radians(lat)))))))
    inds=tree.query(box(lon-rlon,lat-rlat,lon+rlon,lat+rlat))
    if not len(inds):
        return None
    candidates=[]
    for idx in inds:
        poly=geoms[int(idx)]
        if poly.is_empty:continue
        nearest=nearest_points(point,poly)[1]
        _,_,metres=GEOD.inv(lon,lat,nearest.x,nearest.y)
        if metres<=5500:
            candidates.append((metres,str(land_ids[int(idx)]),int(idx)))
    if not candidates:return None
    candidates.sort(key=lambda q:(q[0],int(q[1]),q[2]))
    dist,id,index=candidates[0]
    if dist>5000:
        return None
    # If two original islands/coasts are equally close, the port assignment
    # may be inherently ambiguous (harbor in channel between landmasses).
    alternatives={x[1] for x in candidates if x[0]<=dist+250}
    return {"nearest_original_sibling_id":id,"distance_to_GSHHG_land_km":dist/1000,
            "candidate_landmasses_within_250m_of_nearest":len(alternatives),
            "ambiguous_nearest_landmass":len(alternatives)>1}

def compute(archive,work):
    rows,eligible,frame_sha=read_frame()
    geodata,digest=source_shape(archive,work)
    ids=geodata["landmass_id"].tolist()
    geom=np.asarray(geodata.geometry.to_numpy(),dtype=object)
    tree=STRtree(geom)
    wpi,original_rows=all_wpi()
    output=[];counts=Counter()
    for port in wpi:
        match=nearest_original_landmass(port,tree,geom,ids)
        row={**port}
        row.update(match or {"nearest_original_sibling_id":"",
            "distance_to_GSHHG_land_km":"",
            "candidate_landmasses_within_250m_of_nearest":0,
            "ambiguous_nearest_landmass":False})
        island=bool(match and match["nearest_original_sibling_id"] in eligible)
        row["assigned_candidate_GSHHG_island_area_ge20"]=str(island)
        row["catalogued_2019_WPI_not_historical_ship_density"]=True
        output.append(row)
        counts[("assigned_island" if island else "not_in_ge20_island_frame")]+=1
        if match and match["ambiguous_nearest_landmass"]:
            counts["ambiguous_nearest"]+=1
    all_island=defaultdict(list)
    for x in output:
        if x["assigned_candidate_GSHHG_island_area_ge20"]=="True":
            all_island[x["nearest_original_sibling_id"]].append(x)
    island_rows=[]
    for island in rows:
        ports=all_island[island["landmass_id"]]
        island_rows.append({
            "landmass_id":island["landmass_id"],
            "area_km2":island["area_km2_equal_area"],
            "WPI_2019_assigned_port_count":len(ports),
            "WPI_2019_port_id_list":";".join(str(x["port_id"]) for x in ports),
            "WPI_2019_nearshore_assigned_any":str(bool(ports)),
            "original_17_port_cohort":island["in_original_17_historical_port_botanical_cohort"],
            "WPI_not_shipping_volume":True,
        })
    assert len(island_rows)==len(rows)
    assert sum(x["WPI_2019_assigned_port_count"] for x in island_rows)==counts["assigned_island"]
    HERE.mkdir(exist_ok=True,parents=True)
    def save_csv(name,table):
        with (HERE/name).open("w",encoding="utf8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(table[0]))
            w.writeheader();w.writerows(table)
    save_csv("global_WPI_ports_assigned_to_GSHHG.csv",output)
    save_csv("global_island_port_opportunity.csv",island_rows)
    n_with_ports=sum(x["WPI_2019_assigned_port_count"]>0 for x in island_rows)
    summary={
        "schema":"chocho_global_2019_WPI_island_sibling_assignment_v1",
        "status":"GLOBAL_PORT_PROXIMITY_LOCATOR_NOT_SHIPPING_OR_BUTTERFLY_RANGE_CAUSAL_EVIDENCE",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "source_GSHHG_SHA256":digest,"source_GSHHG_geographic_frame_blob":frame_sha,
        "source_WPI_query":WPI,"WPI_feature_rows":original_rows,"WPI_valid_coordinates":len(wpi),
        "GSHHG_area_threshold_min_km2":AREA_MIN,"GSHHG_continental_exclusion_km2":AREA_MAX,
        "global_island_candidate_count":len(island_rows),
        "GSHHG_islands_with_at_least_one_assigned_modern_WPI_port":n_with_ports,
        "GSHHG_islands_without_indexed_WPI_ports":len(island_rows)-n_with_ports,
        "WPI_port_count_assigned_to_target_islands":counts["assigned_island"],
        "WPI_ports_with_no_ge20_GSHHG_island_assignment":counts["not_in_ge20_island_frame"],
        "WPI_port_identity_ambiguous_within_250m":counts["ambiguous_nearest"],
        "candidate_butterfly_species_for_expansion":225,
        "possible_species_island_at_risk_units":225*len(island_rows),
        "limits":[
            "WPI geography is from 2019 and does not quantify 2001–2010 ship cargo, port opening or introduction timing.",
            "Nearest source coastline within 5km may attach an offshore harbor to the wrong island in narrow channels; ambiguous source counts are reported.",
            "GSHHG high-resolution L1 landmass IDs are not equivalent to established mainland-distance barriers, ocean crossings or true remote island categories.",
            "WPI-unlisted does not mean an island lacked working historical ports or received no ships.",
            "Only butterflies' time-indexed observations, plant species-specific historic host supply and matched controls can address range expansion hypothesis."
        ]
    }
    (HERE/"global_WPI_island_port_receipt.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:summary[k] for k in (
        "WPI_valid_coordinates","global_island_candidate_count",
        "GSHHG_islands_with_at_least_one_assigned_modern_WPI_port",
        "GSHHG_islands_without_indexed_WPI_ports",
        "WPI_port_count_assigned_to_target_islands",
        "possible_species_island_at_risk_units")},ensure_ascii=False,indent=2),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    args=p.parse_args()
    compute(args.zip,args.work)

if __name__=="__main__":main()

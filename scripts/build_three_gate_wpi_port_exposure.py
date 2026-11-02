#!/usr/bin/env python3
"""Q4 port-access proxy: published WPI points to frozen WGSRPD3 island regions.

Modern snapshot counts ports, NOT historic ship traffic or entry paths.
Require identity validation and unique geographic assignment. Keep close-coast
ties as ambiguous instead of attributing ports to neighboring islands.
"""
from __future__ import annotations
import argparse,csv,json,math,time,urllib.parse,urllib.request
from collections import Counter,defaultdict
from pathlib import Path
from shapely.geometry import Point,shape
from shapely.strtree import STRtree
from shapely.ops import transform
from pyproj import Transformer

WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
FIELDS="OBJECTID,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE,INDEX_NO"
AGENT="chocho-three-gates-port-input/0.1"

def request(params):
    url=WPI+"?"+urllib.parse.urlencode(params)
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"}),timeout=90) as r:
        d=json.load(r)
    if "error" in d:raise RuntimeError("WPI error: "+str(d["error"])[:400])
    return d

def readcsv(path):
    with path.open(newline="",encoding="utf-8") as h:yield from csv.DictReader(h)

def save(path,rows):
    with path.open("w",newline="",encoding="utf8") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--island-metrics",type=Path,required=True)
    p.add_argument("--wgsrpd3-geojson",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    allareas={r["area_code_l3"]:r for r in readcsv(a.island_metrics)}
    if len(allareas)!=369:raise RuntimeError("WGSRPD Q1 frozen region census has changed")
    features=json.loads(a.wgsrpd3_geojson.read_text(encoding="utf-8"))["features"]
    shapes=[]
    ids=[]
    for ft in features:
        key=str((ft.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if key not in allareas or not ft.get("geometry"):continue
        geom=shape(ft["geometry"])
        if geom.is_empty:continue
        shapes.append(geom)
        ids.append(key)
    if set(ids)!=set(allareas):raise RuntimeError("Not all 369 WGSRPD polygons available")
    tree=STRtree(shapes)
    xform=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform
    projected=[transform(xform,g) for g in shapes]
    p_tree=STRtree(projected)
    ids_raw=request({"where":"1=1","returnIdsOnly":"true","f":"json"})
    objids=sorted(set(map(int,ids_raw.get("objectIds") or [])))
    if len(objids)<2000:raise RuntimeError("WPI mirror unexpectedly small")
    ports=[]
    for start in range(0,len(objids),200):
        selected=objids[start:start+200]
        page=request({"objectIds":",".join(map(str,selected)),"outFields":FIELDS,
                      "returnGeometry":"false","f":"json"})
        for feature in page.get("features") or []:
            r=feature["attributes"]
            oid=r.get("OBJECTID")
            if oid is None:continue
            ports.append(r)
        time.sleep(.15)
    if len(ports)!=len(objids) or len({int(r["OBJECTID"]) for r in ports})!=len(objids):
        raise RuntimeError(f"WPI incomplete/duplicated: {len(ports)} vs {len(objids)}")
    matched=[]
    diagnostics=Counter()
    by_area=Counter()
    for row in ports:
        oid=int(row["OBJECTID"])
        try:
            lat,lon=float(row["LATITUDE"]),float(row["LONGITUDE"])
            if not -90<=lat<=90 or not -180<=lon<=180:raise ValueError("invalid coordinate")
        except (TypeError,ValueError):
            diagnostics["bad_coordinates"]+=1;continue
        pt=Point(lon,lat)
        candidates={ids[int(i)] for i in tree.query(pt) if shapes[int(i)].covers(pt)}
        method="direct_polygon"
        if len(candidates)!=1:
            ppt=transform(xform,pt)
            neighbors=[int(i) for i in p_tree.query(ppt.buffer(5000))
                       if projected[int(i)].distance(ppt)<=5000]
            nearest=sorted([(projected[i].distance(ppt),ids[i]) for i in neighbors])
            if nearest:
                m=nearest[0][0]
                candidates={name for d,name in nearest if d-m<1000}
                method="nearest_5km_coastal"
        if len(candidates)!=1:
            diagnostics["unmapped_or_ambiguous"]+=1
            continue
        code=next(iter(candidates))
        by_area[code]+=1
        diagnostics[method]+=1
        matched.append({"wpi_objectid":oid,"wgsrpd3_code":code,
                        "geographic_group":allareas[code]["class"],
                        "port_name":row.get("PORT_NAME") or "",
                        "index_no":row.get("INDEX_NO") or "",
                        "assignment_method":method})
    joined=[]
    for code,r in sorted(allareas.items()):
        joined.append({
            "area_code_l3":code,
            "region_name":r["region_name"],
            "level1":r["level1"],
            "geographic_group":r["class"],
            "WPI_port_count":by_area[code],
            "has_indexed_port":int(by_area[code]>0),
            "source_WPI_mirror_last_edited":"2023-08-09",
        })
    save(a.out_dir/"wgsrpd3_port_counts.csv",joined)
    save(a.out_dir/"wpi_port_to_wgsrpd3_assignments.csv",matched)
    islands=[r for r in joined if r["geographic_group"]=="island"]
    result={
      "schema":"chocho_three_gate_wpi_port_island_exposure_v0.1",
      "source":"ArcGIS public WPI point mirror, not an actual historic ship movement ledger",
      "source_url":"https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0",
      "source_mirror_last_edited_utc":"2023-08-09",
      "all_public_WPI_ports":len(ports),
      "matched_port_positions":len(matched),
      "location_diagnostics":dict(diagnostics),
      "frozen_WGSRPD3_regions":len(joined),
      "frozen_strict_island_regions":len(islands),
      "strict_islands_with_at_least_one_WPI_port":sum(x["WPI_port_count"]>0 for x in islands),
      "strict_island_ports_total":sum(x["WPI_port_count"] for x in islands),
      "top_island_port_regions":sorted(
        [{"code":r["area_code_l3"],"name":r["region_name"],"n":r["WPI_port_count"]}
         for r in islands],key=lambda x:-x["n"])[:12],
      "limits":[
       "Modern published WPI indexed port count is not ship visitation, supply of larvae on plants, or port activity at biological first arrival year.",
       "Port number is strongly related to area, urbanization, sampling intensity and historical geography; these require controls.",
       "Coastal ports within 5 km of WGSRPD3 polygon boundaries are assigned only if unambiguous within a 1-km nearest-margin rule.",
       "WPI mirror last edit 2023 differs from original source-locked WPI 2019 used in earlier island-port analysis and cannot be directly compared to it.",
       "World Geographic Scheme regions often aggregate several physical islands; ports do not resolve intra-archipelago biological histories.",
       "Wing×port interaction requires species×region observation risk set and species-specific detectability, not available from this port layer alone."
      ]
    }
    (a.out_dir/"wpi_port_exposure_audit.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("Q4_PORT_INPUT",json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Recover *actual* historic WPI point-source coordinates from NOAA ArcGIS.

The same NOAA nccos/PIRO_DigitalAtlas MapServer/30 WPI2019 endpoint and
OBJECTID used by original chocho source are used; no inferred geography.
Any changed count, duplicate ID or unresolved coordinates fails closed.
Save current source receipt; equality with 2026 archived 51-island distances
MUST still be independently validated before use.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse,urllib.request
from pathlib import Path

BASE="https://gis.ngdc.noaa.gov/arcgis/rest/services/nccos/PIRO_DigitalAtlas/MapServer/30"
EXPECTED_WPI_COUNT=3630

def json_get(path,params):
    u=BASE+path+"?"+urllib.parse.urlencode(params)
    err=None
    for n in range(4):
        try:
            request=urllib.request.Request(u,headers={
                "Accept":"application/json","User-Agent":"chocho-island-2019-WPI-source-audit/1.0"})
            with urllib.request.urlopen(request,timeout=60) as f:
                body=f.read(18000000)
            x=json.loads(body.decode("utf-8-sig"))
            if "error" in x:raise ValueError(str(x["error"])[:300])
            return x,hashlib.sha256(body).hexdigest()
        except Exception as ex:
            err=repr(ex);time.sleep(2*(n+1))
    raise RuntimeError("Official NOAA 2019 WPI endpoint unavailable: "+str(err))

def run(geojson_path,receipt_path):
    info,sha=json_get("",{"f":"json"})
    src={"source_url":BASE,"source_original_publication":"NGA World Port Index 2019 hosted in NOAA PIRO DigitalAtlas",
         "metadata_sha256":sha,
         "metadata_title":info.get("name"),"geometry_type":info.get("geometryType"),
         "max_record_count":info.get("maxRecordCount"),
         "spatial_reference":info.get("extent",{}).get("spatialReference")}
    if str(info.get("name")).lower()!="ports" or info.get("geometryType")!="esriGeometryPoint":
        raise ValueError("NOAA service is not the expected WPI Ports point feature layer")
    count,sha_count=json_get("/query",{"where":"1=1","returnCountOnly":"true","f":"json"})
    total=int(count.get("count",0))
    src["live_port_record_count"]=total
    src["record_count_response_sha256"]=sha_count
    if total!=EXPECTED_WPI_COUNT:
        src["source_changed_gating_failure"]=f"Official 2019 WPI had {EXPECTED_WPI_COUNT}, live {total}"
        receipt_path.parent.mkdir(parents=True,exist_ok=True)
        receipt_path.write_text(json.dumps(src,indent=2)+"\n")
        raise ValueError(src["source_changed_gating_failure"])
    all_features=[]; pages=[]
    for start in range(0,total,900):
        payload,hash_=json_get("/query",{
            "where":"1=1","outFields":"*","outSR":"4326",
            "returnGeometry":"true","orderByFields":"OBJECTID ASC",
            "resultOffset":str(start),"resultRecordCount":"900","f":"json"})
        rows=payload.get("features",[])
        if len(rows)!=min(900,total-start):
            raise ValueError(f"WPI page incomplete at {start} got {len(rows)}")
        for row in rows:
            xy=row.get("geometry") or {}
            attrs=row.get("attributes") or {}
            if "OBJECTID" not in attrs or "x" not in xy or "y" not in xy:
                raise ValueError("Original WPI coordinate record lacks OBJECTID and WGS84 x/y")
            i=int(attrs["OBJECTID"]);x=float(xy["x"]);y=float(xy["y"])
            if not (-180<=x<=180 and -90<=y<=90):
                raise ValueError("Invalid port lon/lat; outSR=4326 not honored")
            all_features.append({"type":"Feature","properties":attrs,
                                  "geometry":{"type":"Point","coordinates":[x,y]}})
        pages.append({"offset":start,"count":len(rows),"sha256":hash_})
    ids=[int(x["properties"]["OBJECTID"]) for x in all_features]
    if len(set(ids))!=total or len(ids)!=total:
        raise ValueError("Nonunique or missing WPI objects")
    all_features.sort(key=lambda v:int(v["properties"]["OBJECTID"]))
    obj={"type":"FeatureCollection","name":"NGA_WPI2019_NOAA_Ports30",
         "features":all_features}
    geojson_path.parent.mkdir(parents=True,exist_ok=True)
    payload=json.dumps(obj,ensure_ascii=False,separators=(",",":"))+"\n"
    geojson_path.write_text(payload,encoding="utf-8")
    src.update({"source_recovered_points":len(all_features),
      "port_id_minmax":[min(ids),max(ids)],
      "georeferenced_geojson_sha256":hashlib.sha256(payload.encode()).hexdigest(),
      "feature_pages":pages,
      "pending_51_landmass_original_distance_validation":True,
      "source_version_cannot_be_assumed_historically_unchanged":True})
    receipt_path.parent.mkdir(parents=True,exist_ok=True)
    receipt_path.write_text(json.dumps(src,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ORIGINAL_NOAA_WPI_2019_PORT_LOCATIONS_REAL_SOURCE",json.dumps(
        {k:v for k,v in src.items() if k!="feature_pages"},sort_keys=True),flush=True)
    return src

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-geojson",type=Path,required=True)
    p.add_argument("--out-receipt",type=Path,required=True)
    a=p.parse_args()
    run(a.out_geojson,a.out_receipt)
if __name__=="__main__":main()

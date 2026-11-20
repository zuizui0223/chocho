#!/usr/bin/env python3
"""Compare all 1596 complete source WPI 2019 point positions with public mirror.

Only publish a source SHA256 check; never publish the original private WPI
GeoJSON or silently change historical geographies. Source 1596 coordinates
were normalized by integer OBJECTID and 7-decimal lat/lon, sorted by OBJECTID.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import URLError,HTTPError

ORIGINAL="3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"
PORTS=Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
PORT_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"
PRIVATE_ORIGINAL_1596_GEO_SHA256="3c10fa14a18f1375ca2c90e0aecd91eb5c0973b598c1c9ec4e9333f5d8092004"
ARCGIS="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"

def get_json(params):
    url=ARCGIS+"?"+urlencode(params)
    last=None
    for retry in range(5):
        try:
            req=Request(url,headers={"User-Agent":"chocho-WPI-2019-original1596-source-20261120",
                      "Accept":"application/json"})
            with urlopen(req,timeout=60) as response:
                if response.status!=200:raise RuntimeError("Public WPI mirror HTTP non-200")
                j=json.load(response)
            if not isinstance(j,dict) or j.get("error"):
                raise ValueError("ArcGIS WPI source JSON error: "+str(j.get("error"))[:150])
            if j.get("exceededTransferLimit"):
                raise RuntimeError("ArcGIS WPI mirror response truncated")
            return j
        except (TimeoutError,HTTPError,URLError,ValueError,RuntimeError) as e:
            last=e
            if isinstance(e,HTTPError) and e.code not in (408,429,500,502,503,504):break
            if retry<4:time.sleep(min(12,2**retry))
    raise RuntimeError(f"WPI source unavailable, not zero records: {type(last).__name__}: {str(last)[:170]}")

def run(root:Path):
    data=(root/PORTS).read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if blob!=PORT_BLOB:raise ValueError("Original 3630 source WPI port table changed")
    rows=list(csv.DictReader(data.decode("utf8").splitlines()))
    ids=sorted(int(r["port_id"]) for r in rows if r["query_complete"]=="True")
    if len(rows)!=3630 or len(ids)!=1596 or len(set(ids))!=1596:
        raise ValueError("Original WPI source complete-port cohort changed")
    coords={}
    for ix,start in enumerate(range(0,len(ids),65),start=1):
        block=ids[start:start+65]
        j=get_json({"where":"OBJECTID IN ("+",".join(map(str,block))+")",
                    "outFields":"OBJECTID,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
                    "returnGeometry":"false","f":"json"})
        records=j.get("features",[])
        if len(records)!=len(block):raise ValueError(f"Missing WPI mirror objects batch {ix}: {len(records)} != {len(block)}")
        for f in records:
            a=f["attributes"]
            k=int(a["OBJECTID"])
            if k not in block or k in coords:raise ValueError("Duplicate/foreign mirrored WPI source ID")
            lat,lon=float(a["LATITUDE"]),float(a["LONGITUDE"])
            if not (-90<=lat<=90 and -180<=lon<=180):raise ValueError("Invalid mirrored coordinates")
            coords[k]=(lat,lon)
        print("WPI1596_COORDINATE_SOURCE_BATCH",ix,"n",len(block),flush=True)
        time.sleep(.12)
    if set(coords)!=set(ids):raise ValueError("Source WPI mirror is not a complete cohort")
    raw="port_id,latitude_7dp,longitude_7dp\n"+"".join(
      f"{k},{coords[k][0]:.7f},{coords[k][1]:.7f}\n" for k in ids)
    digest=hashlib.sha256(raw.encode("utf8")).hexdigest()
    if digest!=PRIVATE_ORIGINAL_1596_GEO_SHA256:
        raise ValueError("Public WPI mirror diverged from private original 2019 point SHA256; no full port requery authorized. Mirror SHA256 "+digest)
    return {"schema":"chocho_original1596_WPI2019_private_source_vs_public_arcgis_coordinate_sha_v1",
       "status":"ALL_ORIGINAL_1596_PORT_ID_COORDINATES_REPRODUCED_FROM_PUBLIC_WPI_MIRROR",
       "source_original_commit":ORIGINAL,
       "full_original_port_candidates":len(rows),
       "original_complete_plant_query_ports":len(ids),
       "matched_WPI_OBJECTID_and_coords":len(coords),
       "public_WPI_endpoint":ARCGIS,
       "original_private_WPI_1596_7dp_coordinate_SHA256":PRIVATE_ORIGINAL_1596_GEO_SHA256,
       "current_public_WPI_1596_7dp_coordinate_SHA256":digest,
       "original_private_WPI_geojson_NOT_UPLOADED":True,
       "GBIF_individual_occurrence_provenance_not_recovered":True,
       "only_port_geography_not_original_plant_query_buffer_hash":True,
       "checked_at_utc":datetime.now(timezone.utc).isoformat()}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    receipt=run(args.source_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("WPI_ORIGINAL1596_SOURCE_GEO_MIRROR_CHECK_PASS",
          receipt["current_public_WPI_1596_7dp_coordinate_SHA256"])
if __name__=="__main__":
    main()

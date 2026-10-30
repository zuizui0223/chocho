#!/usr/bin/env python3
"""Audit ALL 2025 UN-OCHA/NGA WPI ports using original numeric lat/lon attributes.

2025 port coordinates are NOT ship calls, historic transport volumes, port
opening times or butterflies arriving. No raw port or exact coords exported.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time
from pathlib import Path
from collections import Counter
from urllib.parse import urlencode
from audit_island_port_original_sources_20261029 import WPI,request

SCHEMA="chocho_WPI202511_original_3689_port_attribute_geography_v01"
EXPECTED=3689
ORIGINAL_METADATA_SHA="a2729ba58fcf5a74af65fe9f3460bfe3d1c5d5ebd43b05c383bfd52d59369a42"
ORIGINAL_IDS_SHA="8aae61cbba187c5a233188c6aa3a3f905e5b6b64293236e56d111dcb30bfb6d8"

def validate(rows,ids):
    expected=set(ids)
    got=set()
    valid=0;invalid=0;coasts=Counter();quarter=Counter();duplicates=Counter()
    source_numeric_zero=0
    for row in rows:
        a=row.get("attributes")
        if not isinstance(a,dict):raise ValueError("Original WPI attributes absent")
        k={str(x).lower():y for x,y in a.items()}
        ident=k.get("objectid")
        if ident not in expected or ident in got:
            raise ValueError("Invalid/duplicated WPI original port object id")
        got.add(ident)
        try:
            lat=float(k["latitude"]);lon=float(k["longitude"])
        except (ValueError,TypeError,KeyError):
            invalid+=1;continue
        if not(math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180):
            invalid+=1;continue
        valid+=1
        if lat==0 and lon==0:source_numeric_zero+=1
        band=(math.floor(lat/30)*30,math.floor(lon/60)*60)
        quarter[f"{band[0]}:{band[1]}"]+=1
        country=k.get("country")
        if isinstance(country,str) and len(country.strip())<=120:
            coasts[country.strip()]+=1
        # Coordinate duplicates may represent close but distinct registered ports;
        # count and report them, never infer actual shipment redundancy.
        duplicates[(round(lat,4),round(lon,4))]+=1
    if got!=expected:
        raise ValueError(f"Original WPI 2025 incomplete port object IDs {len(got)}/{len(expected)}")
    if invalid>max(25,math.ceil(.02*len(ids))):
        raise ValueError(f"WPI original-coordinate invalid fraction exceeds 2%: {invalid}/{len(ids)}")
    return {
      "source_object_ids_complete":len(got),
      "valid_original_WGS84_attribute_coordinates":valid,
      "missing_or_invalid_attribute_coordinates":invalid,
      "explicit_zerozero_coordinates":source_numeric_zero,
      "coordinate_cells_with_duplicate_four_decimal_locations":sum(v>1 for v in duplicates.values()),
      "WPI_reported_country_name_categories":len(coasts),
      "hemisphere_coarse_30x60_degree_band_counts":dict(sorted(quarter.items())),
      "raw_port_names_and_coordinates_or_IDs_exported":False,
      "original_source_2025_ports_are_NOT_trade_or_ship_arrival_counts":True,
    }

def run(output):
    meta=request(WPI+"?f=json")
    rawids=request(WPI+"/query?"+urlencode({"f":"json","where":"1=1","returnIdsOnly":"true"}))
    if hashlib.sha256(meta).hexdigest()!=ORIGINAL_METADATA_SHA or hashlib.sha256(rawids).hexdigest()!=ORIGINAL_IDS_SHA:
        raise ValueError("Original 2025 WPI source metadata or ID roster drifted: do not replace silently")
    m=json.loads(meta);i=json.loads(rawids)
    if m.get("geometryType")!="esriGeometryPoint":raise ValueError("Original port source not point")
    original=i.get("objectIds",[])
    if not isinstance(original,list) or len(original)!=EXPECTED or len(set(original))!=EXPECTED:
        raise ValueError(f"Original source roster 3689 changed, observed {len(original)}")
    requested=sorted(original)
    rows=[]
    digests=[]
    for offset in range(0,len(requested),200):
        chunk=requested[offset:offset+200]
        url=WPI+"/query?"+urlencode({"f":"json","where":"1=1","objectIds":",".join(map(str,chunk)),
            "returnGeometry":"false","outFields":"objectid,latitude,longitude,country,index_no"})
        response=request(url,1_500_000)
        digest=hashlib.sha256(response).hexdigest()
        j=json.loads(response)
        if "error" in j or not isinstance(j.get("features"),list):
            raise ValueError(f"Original WPI API rejected paginated source chunk {offset//200}")
        if len(j["features"])!=len(chunk):
            raise ValueError(f"Original WPI port chunk incomplete {offset}//200")
        rows.extend(j["features"]);digests.append(digest)
        if (offset//200+1)%5==0:print(json.dumps({"WPI_source_original_port_rows_acquired":len(rows)}),flush=True)
        time.sleep(.1)
    stats=validate(rows,original)
    out={
      "schema":SCHEMA,"status":"WPI_2025_COMPLETE_ORIGINAL_3689_PORT_ATTRIBUTE_COORDINATES_VERIFIED",
      "source_uri":WPI,"original_WPI_metadata_sha256":ORIGINAL_METADATA_SHA,
      "original_WPI_objectID_source_body_sha256":ORIGINAL_IDS_SHA,
      "verified_source_pagination_chunks":len(digests),
      "SHA256_of_canonical_sorted_original_page_SHA256s":
          hashlib.sha256("\n".join(sorted(digests)).encode()).hexdigest(),
      "all_port_2025_original_data":stats,
      "no_time_stamped_maritime_traffic_or_shipping_volume_model":True,
      "no_butterfly_establishment_or_isolation_causal_result":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["status"],"valid":stats["valid_original_WGS84_attribute_coordinates"],
        "invalid":stats["missing_or_invalid_attribute_coordinates"],
        "duplicate_location_cells":stats["coordinate_cells_with_duplicate_four_decimal_locations"]},sort_keys=True))
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    run(p.parse_args().output)

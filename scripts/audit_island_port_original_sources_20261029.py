#!/usr/bin/env python3
"""Read-only ORIGINAL Natural Earth 10m land and OCHA/NGA 2025 WPI source probe.

This probes stable geographic SOURCE coverage only. Source port points,
historical route/island names, detailed private observed butterfly coordinates,
H3-host support and raw features are NOT exported. Exact SHA and source
schema/safe count are reported, without claiming port-traffic exposure.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen

NE="https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip"
WPI="https://gis.unocha.org/server/rest/services/Hosted/global_world_seaport_index_202511/FeatureServer/0"
SCHEMA="chocho_island_port_NaturalEarth10m_UNOCHA_WPI202511_original_source_gate_v01"
MAX_META_BYTES=3_000_000

def request(url,max_bytes=MAX_META_BYTES):
    req=Request(url,headers={"User-Agent":"chocho-geographic-source-audit/1.0",
                              "Accept":"application/json,application/zip,application/octet-stream,*/*"})
    with urlopen(req,timeout=75) as r:
        raw=r.read(max_bytes+1)
        if r.status!=200:raise ValueError("Source HTTP status not 200")
    if len(raw)>max_bytes:raise ValueError("Source exceeds predeclared safety bound")
    return raw

def validate_ne_land(source):
    if not source.startswith(b"PK\x03\x04"):raise ValueError("NE original must be ZIP")
    with zipfile.ZipFile(io.BytesIO(source)) as z:
        if z.testzip() is not None:raise ValueError("Natural Earth ZIP corrupt")
        names=z.namelist()
        for suffix in (".shp",".shx",".dbf",".prj"):
            if len([n for n in names if n.endswith("ne_10m_land"+suffix)])!=1:
                raise ValueError("Incomplete Natural Earth land shapefile source")
        if sum(x.file_size for x in z.infolist())>40_000_000:
            raise ValueError("Source land ZIP uncompressed too large")
    return {"source":"Natural Earth 10m land polygons (not exact island objects)",
            "source_uri":NE,"original_zip_bytes":len(source),
            "source_original_zip_sha256":hashlib.sha256(source).hexdigest(),
            "shapefile_parts_verified":["shp","shx","dbf","prj"],
            "island_connected_components_not_yet_validated":True}

def check_wpi_metadata(metadata,ids,sample):
    if "error" in metadata:raise ValueError("WPI layer source denied")
    if metadata.get("type") not in ("Feature Layer","Feature Service","FeatureServer") and metadata.get("geometryType")!="esriGeometryPoint":
        raise ValueError("WPI original layer not point geography")
    fields=[f.get("name","") for f in metadata.get("fields",[]) if isinstance(f,dict)]
    if len(fields)<5:raise ValueError("WPI point schema insufficient")
    if "error" in ids or "error" in sample:raise ValueError("WPI data query inaccessible")
    key=ids.get("objectIdFieldName")
    keys=ids.get("objectIds")
    if not isinstance(keys,list) or not 1000<len(keys)<10000:
        raise ValueError("WPI original port index source count missing/outside expected range")
    if not isinstance(sample.get("features"),list) or not sample["features"]:
        raise ValueError("WPI sample data lacks original locations")
    def pair_ok(lat,lon):
        try:
            a=float(lat);b=float(lon)
            return -90<=a<=90 and -180<=b<=180 and a==a and b==b
        except (TypeError,ValueError):
            return False
    geometry_ok=0
    attribute_ok=0
    attribute_conflict=0
    sample_n=len(sample["features"])
    for row in sample["features"]:
        if not isinstance(row,dict):continue
        geom=row.get("geometry") or {}
        attr={str(k).upper():v for k,v in (row.get("attributes") or {}).items()}
        geom_valid=pair_ok(geom.get("y"),geom.get("x"))
        attr_valid=pair_ok(attr.get("LATITUDE"),attr.get("LONGITUDE"))
        geometry_ok+=int(geom_valid)
        attribute_ok+=int(attr_valid)
        if geom_valid and attr_valid:
            attribute_conflict+=int(abs(float(geom["y"])-float(attr["LATITUDE"]))>1. or
                                    abs(float(geom["x"])-float(attr["LONGITUDE"]))>1.)
    if max(geometry_ok,attribute_ok)<min(5,sample_n):
        raise ValueError("WPI both spatial geometry and original LATITUDE/LONGITUDE source attributes invalid")
    # In this global projected layer the FeatureServer's outSR=4326 may not
    # preserve all global points correctly. Original WPI decimal coordinates
    # are audited separately; refuse to silently replace them with dubious
    # transformed geometry.
    if attribute_ok>=min(5,sample_n):
        policy="ORIGINAL_WPI_LATITUDE_LONGITUDE_ATTRIBUTES"
    else:
        policy="ARCGIS_OUTSR_4326_POINTS"
    # Source conflicts remain visible and must be resolved before *all* points
    # are used as geodetic sites in a downstream port-distance calculation.
    return {"port_service_uri":WPI,
            "layer_name":metadata.get("name",""),
            "layer_description_n_source_fields":len(fields),
            "source_feature_count":len(keys),
            "original_objectid_field":key,
            "source_fields_available":sorted(fields),
            "sample_points":sample_n,
            "sample_valid_ArcGIS_4326_geometry":geometry_ok,
            "sample_valid_original_attribute_latlon":attribute_ok,
            "sample_geometry_vs_attribute_mismatch_over_1deg":attribute_conflict,
            "preferred_position_source":policy,
            "last_update_202511_is_not_opening_date":True,
            "real_commodity_flow_or_ship_calls_not_recorded":True}

def run(output,manifest_path):
    h3=json.loads(manifest_path.read_text(encoding="utf-8"))
    if h3.get("schema")!="chocho_global_h3r4_WGSRPD3_fractional_geometry_bridge_20261029_v01" or (
       h3["coverage"].get("cells_any_tdwg_land_intersection")!=89161):
        raise ValueError("PR100 source land crosswalk reference changed")
    result={"schema":SCHEMA,
            "source_h3_geometry_bridge_sha256":hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            "h3_source_original_land_cells":89161,
            "h3_no_centroid_nearshore_cells":h3["coverage"]["coastal_cells_without_tdwg_centroid"],
            "original_source_data_or_port_coordinates_serialized":False,
            "does_not_impute_island_discontinuity_from_h3_center":True,
            "no_butterfly_island_port_causal_analysis":True}
    try:
        natural=request(NE,10_000_000)
        result["land"]=validate_ne_land(natural)
    except (OSError,ValueError,zipfile.BadZipFile) as e:
        result["land_error"]={"type":type(e).__name__,"description":str(e)[:300]}
    try:
        metadata_raw=request(WPI+"?f=json")
        metadata=json.loads(metadata_raw)
        ids_raw=request(WPI+"/query?"+urlencode({"f":"json","where":"1=1","returnIdsOnly":"true"}))
        ids=json.loads(ids_raw)
        oid=ids.get("objectIds",[])
        if not isinstance(oid,list) or len(oid)<2:raise ValueError("WPI returned no source IDs")
        sample_raw=request(WPI+"/query?"+urlencode({"f":"json","where":"1=1",
                         "returnGeometry":"true","outSR":"4326","outFields":"*",
                         "objectIds":",".join(map(str,oid[:10]))}))
        sample=json.loads(sample_raw)
        result["ports"]=check_wpi_metadata(metadata,ids,sample)
        result["ports"]["full_original_coordinates_not_yet_decoded_or_verified"]=True
        result["port_metadata_sha256"]=hashlib.sha256(metadata_raw).hexdigest()
        result["port_ID_set_sha256"]=hashlib.sha256(ids_raw).hexdigest()
    except (OSError,ValueError,KeyError,json.JSONDecodeError) as e:
        result["port_error"]={"type":type(e).__name__,"description":str(e)[:300]}
    result["status"]=(
      "ALL_ORIGINAL_ISLAND_PORT_SOURCES_IDENTIFIED"
      if "land" in result and "ports" in result else "ISLAND_PORT_SOURCE_INCOMPLETE")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],
        "land_digest":result.get("land",{}).get("source_original_zip_sha256"),
        "source_ports":result.get("ports",{}).get("source_feature_count"),
        "port_error":result.get("port_error"),
        "land_error":result.get("land_error")},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    d=run(a.output,a.manifest)
    if d["status"]!="ALL_ORIGINAL_ISLAND_PORT_SOURCES_IDENTIFIED":raise SystemExit(2)

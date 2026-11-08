#!/usr/bin/env python3
"""Bounded 13-port PortWatch registry WPI2019 independent name+country+coordinate audit.

Exactly 13 source-only preselected port IDs appearing as GLSN2018 candidate,
UNCTAD PLSCIQ candidate, but NOT independently WPI-verified in Zenodo2020 AIS source (NOT no actual ship traffic).
At most 2 official public registry queries per WPI centroid, max 15 returned rows,
NEVER any bulk IMF source extraction or maritime edge-list API requests.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,unicodedata,urllib.request,urllib.error
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode

try:
    import pycountry
except ImportError:
    pycountry=None  # optional for paper-only test environments; source-run CI installs it

ROOT=Path(__file__).resolve().parents[1]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
JOIN=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv"
GLSN=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_originalWPI_source_candidate_model_edges.csv"
AIS=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv"
REGISTRY="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0"
IDS=("96","98","146","284","402","408","410","415","416","576","656","720","1516")
HALF_BOX_DEGREES=.35
MAX_PROVIDER_ROWS_PER_BOX=15
MAX_SOURCE_QUERIES=2*len(IDS)
EARTH_KM=6371.0088
MATCH_KM=25.0
PRIOR_VALIDATED={}

def get_source_rows(path):
    with Path(path).open("r",encoding="utf-8",newline="") as file:
        return list(csv.DictReader(file))

def gitblob(path):
    data=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def normalize(source):
    raw=unicodedata.normalize("NFKD",str(source))
    return " ".join(re.findall(r"[A-Z0-9]+",raw.encode("ascii","ignore").decode("ascii").upper()))

def haversine(lat1,lon1,lat2,lon2):
    a,b,c,d=[math.radians(float(s)) for s in (lat1,lon1,lat2,lon2)]
    h=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 2*EARTH_KM*math.asin(math.sqrt(max(0.,min(1.,h))))

def frozen_thirteen():
    if gitblob(WPI)!=WPI_BLOB:
        raise ValueError("Original 194 WPI source pin changed")
    ports=get_source_rows(WPI)
    join=get_source_rows(JOIN)
    if len(ports)!=194 or len(join)!=194 or len({p["island_id"] for p in ports})!=51 or sum(p["included"].lower()=="true" for p in ports)!=182:
        raise ValueError("Original 51 islands or 194 WPI / 182 marine source ports changed")
    port_by_id={str(int(p["port_id"])):p for p in ports}
    if len(port_by_id)!=194:
        raise ValueError("Duplicate original WPI port ID")
    selected={}
    for row in join:
        pid=str(int(row["original_WPI_port_id"]))
        if (row["GLSN2018_status"].startswith("UNIQUE_") and
            row["UNCTAD_PLSCIQ_match_status"].startswith("UNIQUE_") and
            row["Zenodo2020_actual_WPI_port_identity_verified"]!="True"):
            if pid in selected or pid not in port_by_id:
                raise ValueError("Duplicate or unknown cross-source original WPI identity")
            w=port_by_id[pid]
            if (row["island_id"]!=w["island_id"] or row["original_WPI_name"]!=w["port_name"] or
                row["original_WPI_country_ISO2"]!=w["country_code"]):
                raise ValueError("Original source candidate and WPI original port/name disagree")
            selected[pid]=(w,row)
    if set(selected)!=set(IDS):
        raise ValueError("Source-selected 13 GLSN+UNCTAD ports lacking verified 2020 AIS endpoint changed")
    if len({selected[i][0]["island_id"] for i in IDS})!=5:
        raise ValueError("Original 13 source-selected ports do not represent exactly five island groups")
    return selected

def api(params):
    url=REGISTRY+"/query?"+urlencode(params)
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-original13-WPI-source-only-IMF-PortWatch-registry/1.0",
        "Accept":"application/json"})
    receipt={"request_url":url,"http_status":None,"response_bytes":0,"error":None}
    try:
        with urllib.request.urlopen(req,timeout=45) as response:
            raw=response.read(650001)
            receipt["http_status"]=response.status
        if len(raw)>650000:
            raise ValueError("Bounded IMF source registry query returned too many bytes")
        receipt["response_bytes"]=len(raw)
        document=json.loads(raw.decode("utf-8"))
        if not isinstance(document,dict) or document.get("error"):
            raise ValueError("Provider returned ArcGIS error: "+str(document.get("error"))[:200])
        return document,receipt
    except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ConnectionError,ValueError) as exc:
        receipt["error"]={"kind":type(exc).__name__,"detail":str(exc)[:200]}
        return None,receipt

def verify_provider_box(pid,w,expected_iso3,features,count):
    original_lat=float(w["original_latitude"])
    original_lon=float(w["original_longitude"])
    candidate=[]
    full=[]
    for feature in features:
        attrs=(feature or {}).get("attributes") or {}
        row={k:attrs.get(k) for k in ("portid","portname","country","ISO3","lat","lon","LOCODE")}
        try:
            lat=float(row["lat"]);lon=float(row["lon"])
            good_coordinate=math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180
            distance=haversine(original_lat,original_lon,lat,lon) if good_coordinate else None
        except (TypeError,ValueError,OverflowError):
            distance=None
        name=bool(row["portname"]) and normalize(row["portname"])==normalize(w["port_name"])
        country=bool(row["ISO3"]) and row["ISO3"]==expected_iso3
        exact=bool(name and country and distance is not None and distance<=MATCH_KM)
        row.update({"source_name_exact":name,"source_country_ISO3_exact":country,
                    "source_geodesic_km":distance,"all_three_identity_checks":exact})
        full.append(row)
        if exact:
            candidate.append(row)
    complete=(isinstance(count,int) and count<=MAX_PROVIDER_ROWS_PER_BOX and count==len(full))
    if not complete:status="INCOMPLETE_OR_FAILED_REGISTRY_BOX"
    elif len(candidate)==1:status="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH"
    elif len(candidate)>1:status="AMBIGUOUS_MULTIPLE_PORTWATCH_NAMES_IN_BOX"
    else:status="NO_EXACT_SOURCE_PORT_MATCH_IN_COMPLETE_BOX"
    return full,(candidate[0] if status=="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH" else None),status

def write_csv(path,rows,fields):
    with Path(path).open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        writer.writeheader();writer.writerows(rows)

def run(out_dir):
    if pycountry is None:
        raise RuntimeError("pycountry source ISO3 mapping required; dedicated source workflow installs it")
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/s).resolve()) for s in ("data","results","provenance")):
        raise ValueError("Do not overwrite original scientific data")
    source=frozen_thirteen()
    out.mkdir(parents=True,exist_ok=True)
    all_checks=[]
    detailed=[]
    for pid in IDS:
        w,join=source[pid]
        country=pycountry.countries.get(alpha_2=w["country_code"])
        if country is None:
            raise ValueError("No independent ISO2-to-ISO3 mapping for original WPI country")
        lat=float(w["original_latitude"]);lon=float(w["original_longitude"])
        if abs(lat)>89.65 or abs(lon)>179.65:
            raise ValueError("Original WPI source coordinate unsafe for a fixed .35-degree bounded registry box")
        where=(f"lat >= {lat-HALF_BOX_DEGREES:.10f} AND lat <= {lat+HALF_BOX_DEGREES:.10f} "
               f"AND lon >= {lon-HALF_BOX_DEGREES:.10f} AND lon <= {lon+HALF_BOX_DEGREES:.10f}")
        result_count,count_req=api({"where":where,"returnCountOnly":"true","f":"json"})
        n=result_count.get("count") if isinstance(result_count,dict) else None
        if n is not None and not isinstance(n,int):
            raise ValueError("Public registry count-only returned a noninteger")
        sample,sample_req=api({
            "where":where,
            "outFields":"portid,portname,country,ISO3,lat,lon,LOCODE",
            "resultRecordCount":str(MAX_PROVIDER_ROWS_PER_BOX),
            "returnGeometry":"false","f":"json"})
        features=sample.get("features",[]) if isinstance(sample,dict) else []
        rows,candidate,status=verify_provider_box(pid,w,country.alpha_3,features,n)
        if (count_req["error"] or sample_req["error"] or
            (sample or {}).get("exceededTransferLimit") is True):
            status="INCOMPLETE_OR_FAILED_REGISTRY_BOX"
            candidate=None
        if status=="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH" and candidate is None:
            raise ValueError("Confirmed source identity status without a complete unique sample")
        # Retain raw bounded source metadata ONLY in CI artifact; public Git row
        # below includes only independently confirmed WPI metadata or unmatched status.
        detailed.append({"original_WPI_port_id":pid,"where":where,"provider_count":n,
            "count_query":count_req,"bounded_sample_query":sample_req,
            "source_provider_rows":rows,"match_status":status})
        old=PRIOR_VALIDATED.get(pid)
        prior_compatibility=(None if old is None or candidate is None
                             else str(candidate["portid"])==old)
        all_checks.append({
            "original_WPI_port_id":pid,
            "original_island_id":w["island_id"],
            "original_WPI_port_name":w["port_name"],
            "original_country_ISO2":w["country_code"],
            "original_WPI_lat":lat,"original_WPI_lon":lon,
            "GLSN2018_modelled_port_name":join["GLSN2018_provider_port_name"],
            "UNCTAD_historical_port_series":join["UNCTAD_PLSCIQ_match_status"],
            "UNCTAD_2006_2010_usable_quarters":join["UNCTAD_2006_2010_available_quarters"],
            "external_registry_box_status":status,
            "source_bounded_registry_records":len(rows),
            "source_bounded_registry_count":n,
            "verified_PortWatch_portid":candidate["portid"] if candidate else "",
            "verified_PortWatch_portname":candidate["portname"] if candidate else "",
            "verified_PortWatch_ISO3":candidate["ISO3"] if candidate else "",
            "verified_PortWatch_LOCODE":candidate["LOCODE"] if candidate and candidate["LOCODE"] else "",
            "verified_PortWatch_lat":candidate["lat"] if candidate else "",
            "verified_PortWatch_lon":candidate["lon"] if candidate else "",
            "original_WPI_to_PortWatch_distance_km":candidate["source_geodesic_km"] if candidate else "",
            "prior_pilot_provider_id":old or "",
            "prior_pilot_compatible_if_both_verified":prior_compatibility if prior_compatibility is not None else "",
            "original_WPI_independent_geographic_identity_verified":bool(candidate),
            "provider_source_are_not_observed_ship_voyages":True,
        })
        print("PORTWATCH13_SOURCE_IDENTITY_PROGRESS="+json.dumps({
            "id":pid,"name":w["port_name"],"island":w["island_id"],"bounded_count":n,
            "status":status,"provider_portid":candidate["portid"] if candidate else None
        },sort_keys=True),flush=True)
    if len(all_checks)!=13 or len(detailed)!=13:
        raise ValueError("Original frozen 13 source support incomplete")
    source_unique={}
    for r in all_checks:
        portid=r["verified_PortWatch_portid"]
        if portid:
            if portid in source_unique:
                raise ValueError("One source PortWatch portid assigned to two different WPI ports")
            source_unique[portid]=r["original_WPI_port_id"]
        if r["prior_pilot_compatible_if_both_verified"] is False:
            raise ValueError("Four previously verified pilot port source IDs inconsistent; inspect source drift")
    write_csv(out/"WPI13_independent_PortWatch_source_geography_LOCODE.csv",all_checks,list(all_checks[0]))
    (out/"bounded_13_source_query_rows.json").write_text(
        json.dumps(detailed,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    receipt={"schema":"chocho_original51_13_independent_IMF_PortWatch_WPI_port_geoverification_v1",
        "status":"ALL_13_SOURCE_SELECTED_WPI_PORT_BOXES_AUDITED_WITHOUT_GLOBAL_BULK_API",
        "source_registry":REGISTRY,
        "source_provider":"International Monetary Fund PortWatch",
        "source_attribution":"Source: IMF PortWatch; registry source access and usage conditions remain provider governed.",
        "source_license":"IMF public data usage terms; no systematic global download authorized",
        "retrieved_at_utc":datetime.now(timezone.utc).isoformat(),
        "source_WPI_original_git_blob":WPI_BLOB,
        "source_original_WPI_ports":194,
        "source_original_islands":51,
        "source_original_marine_ports":182,
        "source_selected_13_WPI_port_IDs":list(IDS),
        "source_selected_13_landmass_count":len({r["original_island_id"] for r in all_checks}),
        "selected_because_GLSN_UNCTAD_but_Zenodo2020_WPI_endpoint_UNMATCHED":True,
        "maximum_public_provider_queries":MAX_SOURCE_QUERIES,
        "provider_queries_used":len(IDS)*2,
        "max_bounded_provider_features_per_box":MAX_PROVIDER_ROWS_PER_BOX,
        "source_verified_WPI_port_count":len(source_unique),
        "source_verified_original_island_count":len({
            r["original_island_id"] for r in all_checks if r["verified_PortWatch_portid"]}),
        "source_verified_WPI_ports_with_LOCODE":sum(bool(r["verified_PortWatch_LOCODE"]) for r in all_checks),
        "source_unmatched_or_incomplete_source_WPI_ports":sum(not r["verified_PortWatch_portid"] for r in all_checks),
        "source_status_counts":dict(Counter(r["external_registry_box_status"] for r in all_checks)),
        "exact_original_4_prior_pilot_ids_unchanged_when_reverified":all(
            r["prior_pilot_compatible_if_both_verified"] is not False for r in all_checks),
        "no_global_registry_download":True,
        "no_global_maritime_connection_edges_downloaded":True,
        "no_butterfly_fits_or_true_ship_voyage_claim":True,
        "all_13_source_rows_and_actual_query_provenance_SHA256":hashlib.sha256(
          (out/"bounded_13_source_query_rows.json").read_bytes()).hexdigest()}
    (out/"WPI13_PortWatch_source_identity_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PORTWATCH13_INDEPENDENT_GEO_RESULT_JSON="+json.dumps({
        "n":len(IDS),"verified_WPI":len(source_unique),
        "verified_islands":receipt["source_verified_original_island_count"],
        "with_LOCODE":receipt["source_verified_WPI_ports_with_LOCODE"],
        "status":receipt["source_status_counts"],
        "portids":[{"wpi":r["original_WPI_port_id"],"portid":r["verified_PortWatch_portid"],
                    "locode":r["verified_PortWatch_LOCODE"],"distance":r["original_WPI_to_PortWatch_distance_km"]}
                   for r in all_checks if r["verified_PortWatch_portid"]]},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)

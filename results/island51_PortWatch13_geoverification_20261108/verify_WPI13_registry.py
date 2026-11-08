#!/usr/bin/env python3
"""Independent stdlib source-level audit of 13 bounded PortWatch source registry boxes.

Reads full bounded provider rows only within ephemeral CI artifact, not distributed
as an IMF registry corpus. Checks all 13 fixed original WPI sources, each sampled
PortWatch ISO3/name/coordinate identity, provider query provenance, and published
source-only verified WPI/LOCODE ledger. No butterfly outcome values.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse,parse_qs

ROOT=Path(__file__).resolve().parents[2]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
JOIN=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
REGISTRY="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0"
IDS=("96","98","146","284","402","408","410","415","416","576","656","720","1516")
PRIOR_IDS={}
RADIUS=6371.0088
MAX=25.0

def read(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def norm(v):
    x=unicodedata.normalize("NFKD",str(v)).encode("ascii","ignore").decode("ascii").upper()
    return " ".join(re.findall("[A-Z0-9]+",x))

def hav(a,b,c,d):
    lat1,lon1,lat2,lon2=[math.radians(float(v)) for v in (a,b,c,d)]
    h=math.sin((lat2-lat1)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 2*RADIUS*math.asin(math.sqrt(min(1.,max(0.,h))))

def close(a,b,label):
    if not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-8):
        raise ValueError("Original source mismatch: "+label)

def verify(directory):
    dir=Path(directory)
    data=WPI.read_bytes()
    if hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()!=WPI_BLOB:
        raise ValueError("Original WPI2019 source Git blob changed")
    ports=read(WPI)
    prior={str(int(row["port_id"])):row for row in ports}
    j=read(JOIN)
    cross={str(int(row["original_WPI_port_id"])):row for row in j}
    if len(ports)!=194 or len(prior)!=194 or len(j)!=194 or len(cross)!=194:
        raise ValueError("Original 194 source ports missing or duplicated")
    picked={pid for pid,row in cross.items() if row["GLSN2018_status"].startswith("UNIQUE_") and
        row["UNCTAD_PLSCIQ_match_status"].startswith("UNIQUE_") and
        row["Zenodo2020_actual_WPI_port_identity_verified"]!="True"}
    if picked!=set(IDS):
        raise ValueError("13 GLSN/UNCTAD source ports not previously observed by Zenodo2020 AIS selected from altered cohort")
    full=json.loads((dir/"bounded_13_source_query_rows.json").read_text(encoding="utf-8"))
    table=read(dir/"WPI13_independent_PortWatch_source_geography_LOCODE.csv")
    receipt=json.loads((dir/"WPI13_PortWatch_source_identity_receipt.json").read_text(encoding="utf-8"))
    if len(full)!=13 or len(table)!=13 or receipt["provider_queries_used"]!=26 or receipt["source_selected_13_WPI_port_IDs"]!=list(IDS):
        raise ValueError("Original complete bounded 13 sample source evidence missing")
    byid={r["original_WPI_port_id"]:r for r in table}
    if set(byid)!=set(IDS):
        raise ValueError("Archived 13 source port status table ID coverage differs")
    if hashlib.sha256((dir/"bounded_13_source_query_rows.json").read_bytes()).hexdigest()!=receipt["all_13_source_rows_and_actual_query_provenance_SHA256"]:
        raise ValueError("Source sampled registry provider rows/provenance SHA256 mismatch")
    statuses=Counter()
    matched=set()
    matched_islands=set()
    nlocode=0
    all_sample_rows=0
    for i,pid in enumerate(IDS):
        source=full[i]
        published=byid[pid]
        wpi=prior[pid]
        crossrow=cross[pid]
        if source["original_WPI_port_id"]!=pid or published["original_WPI_port_id"]!=pid:
            raise ValueError("Original WPI fixed source order changed")
        if (published["original_island_id"]!=wpi["island_id"] or
            published["original_WPI_port_name"]!=wpi["port_name"] or
            published["original_country_ISO2"]!=wpi["country_code"] or
            published["GLSN2018_modelled_port_name"]!=crossrow["GLSN2018_provider_port_name"] or
            published["UNCTAD_historical_port_series"]!=crossrow["UNCTAD_PLSCIQ_match_status"]):
            raise ValueError("Published 13 source identities differ from original WPI and independent sources")
        close(published["original_WPI_lat"],wpi["original_latitude"],"original 13 WPI latitude")
        close(published["original_WPI_lon"],wpi["original_longitude"],"original 17 WPI longitude")
        provider_count=source["provider_count"]
        rows=source["source_provider_rows"]
        all_sample_rows+=len(rows)
        for role in ("count_query","bounded_sample_query"):
            record=source[role]
            u=urlparse(record["request_url"])
            if u.scheme!="https" or u.netloc!="services9.arcgis.com" or u.path!=urlparse(REGISTRY+"/query").path:
                raise ValueError("PortWatch registry source request host/path changed")
            if len(parse_qs(u.query).get("where",[]))!=1:
                raise ValueError("Original source query omitted bounded WHERE")
        complete=(isinstance(provider_count,int) and provider_count<=15 and provider_count==len(rows)
                  and source["count_query"]["error"] is None and source["bounded_sample_query"]["error"] is None)
        qualifying=[]
        expected_iso3=crossrow["GLSN2018_provider_ISO3"]
        for row in rows:
            if not row.get("portid"):
                raise ValueError("Source sample has missing registry portid")
            name=norm(row.get("portname") or "")==norm(wpi["port_name"])
            country=row.get("ISO3")==expected_iso3
            d=None
            try:
                d=hav(wpi["original_latitude"],wpi["original_longitude"],row["lat"],row["lon"])
            except (TypeError,ValueError,OverflowError):
                pass
            if row["source_name_exact"]!=name or row["source_country_ISO3_exact"]!=country or row["source_geodesic_km"] is None and d is not None:
                raise ValueError("Independent registry source exact name/country/geodesic flag disagrees")
            if d is not None:
                close(row["source_geodesic_km"],d,"provider source port haversine")
            valid=bool(name and country and d is not None and d<=MAX)
            if row["all_three_identity_checks"]!=valid:
                raise ValueError("Geography-only provider source sample was promoted incorrectly")
            if valid:
                qualifying.append(row)
        if not complete: status="INCOMPLETE_OR_FAILED_REGISTRY_BOX"
        elif len(qualifying)==1:status="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH"
        elif len(qualifying)>1:status="AMBIGUOUS_MULTIPLE_PORTWATCH_NAMES_IN_BOX"
        else:status="NO_EXACT_SOURCE_PORT_MATCH_IN_COMPLETE_BOX"
        if status!=source["match_status"] or status!=published["external_registry_box_status"]:
            raise ValueError("Source 17 original provider match/ambiguous/incomplete case disagree")
        if status=="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH":
            original=qualifying[0]
            if (published["verified_PortWatch_portid"]!=str(original["portid"]) or
                published["verified_PortWatch_portname"]!=str(original["portname"]) or
                published["verified_PortWatch_ISO3"]!=str(original["ISO3"]) or
                published["verified_PortWatch_LOCODE"]!=str(original["LOCODE"] or "")):
                raise ValueError("Saved provider portid/LOCODE identity differs")
            close(published["original_WPI_to_PortWatch_distance_km"],original["source_geodesic_km"],"published verified geodesic")
            if published["original_WPI_independent_geographic_identity_verified"].lower()!="true":
                raise ValueError("True source verified port missing verified flag")
            if pid in PRIOR_IDS and original["portid"]!=PRIOR_IDS[pid]:
                raise ValueError("Previous four source-only validated WPI pilot IDs changed")
            matched.add(str(original["portid"]))
            matched_islands.add(wpi["island_id"])
            nlocode+=bool(original["LOCODE"])
        else:
            if (published["verified_PortWatch_portid"] or
                published["original_WPI_independent_geographic_identity_verified"].lower()!="false"):
                raise ValueError("An unverified original source port ID was promoted")
        statuses[status]+=1
        if int(published["source_bounded_registry_records"])!=len(rows):
            raise ValueError("Source public registry row count differs")
        if published["source_bounded_registry_count"] and int(published["source_bounded_registry_count"])!=provider_count:
            raise ValueError("Public PortWatch source count-only result altered")
    if (len(matched)!=receipt["source_verified_WPI_port_count"] or
        len(matched_islands)!=receipt["source_verified_original_island_count"] or
        nlocode!=receipt["source_verified_WPI_ports_with_LOCODE"] or
        statuses!=Counter(receipt["source_status_counts"])):
        raise ValueError("All13 source port matches/island groups/status receipt could not be independently reconstructed")
    if all_sample_rows>13*15 or receipt["no_global_registry_download"] is not True:
        raise ValueError("PortWatch source-only bounded query limit violated")
    result={
        "schema":"chocho_13_IMF_PortWatch_geographic_port_identity_independent_source_audit_v1",
        "status":"ALL13_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED",
        "original_WPI_source_ports":194,
        "source_thirteen_port_selection_verified":13,
        "public_registry_queries_independently_checked":26,
        "max_bounded_provider_example_records":13*15,
        "actual_source_example_records":all_sample_rows,
        "independent_geography_portid_verified":len(matched),
        "verified_original_landmass_groups":len(matched_islands),
        "verified_source_port_LOCODE_present":nlocode,
        "unverified_provider_rows_never_promoted_to_WPI":True,
        "no_bulk_IMF_maritime_network_download":True,
        "no_butterfly_or_seed_cargo_inference":True,
    }
    (dir/"independent_WPI13_source_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("PORTWATCH13_INDEPENDENT_SOURCE_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",required=True,type=Path)
    verify(p.parse_args().results_dir)

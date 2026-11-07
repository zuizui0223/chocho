#!/usr/bin/env python3
"""Probe ONLY the 4x4 PortWatch static maritime source subgraph on preverified WPI ports.

No full network download, no vessel trajectories, no butterfly outcomes.
A PortWatch Spillover Simulator maritime edge is an AIS-derived network
relationship, NOT one confirmed independent observed shipping trip.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

ROOT=Path(__file__).resolve().parents[1]
ROOT_RESULT=ROOT/"results/island51_portwatch_global_connections_20261107"
PROVIDER_SAMPLE=ROOT_RESULT/"WPI5_bounded_registry_source_pilot_receipt.json"
SOURCE_LAYER=("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/"
              "Spillover_Simulator_Maritime_Connections/FeatureServer/0")
WPI_IDS=("28","428","817","1416")
PROVIDER_IDS=("port141","port1153","port647","port254")
MAX_SAMPLE_ROWS=16
REQUEST_FIELDS=("from_id","to_id","from_name","to_name","from_iso3","to_iso3",
                "from_lat","from_lon","to_lat","to_lon","distance","transit_days")

def fixed_samples():
    data=json.loads(PROVIDER_SAMPLE.read_text(encoding="utf-8"))
    if (data.get("status")!="FIVE_ORIGINAL_WPI_SOURCE_PORT_PROVIDER_REGISTRY_SAMPLES_COMPLETE" or
        data.get("original_WPI_source_ports")!=194 or
        data.get("original_source_islands")!=51 or
        data.get("source_WPI_sample_size")!=5):
        raise ValueError("Original WPI5 public provider source pilot missing or altered")
    by_id={x["original_WPI_port_id"]:x for x in data["source_check_rows"]}
    if len(by_id)!=5:raise ValueError("Five WPI identity-check source rows missing")
    selected=[]
    for wpi,portid in zip(WPI_IDS,PROVIDER_IDS):
        p=by_id.get(wpi)
        if (p is None or p["verified_provider_portid"]!=portid or
            p["source_port_identity_verification_status"]!="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH"):
            raise ValueError("Original four WPI source ports no longer have verified PortWatch provider IDs")
        selected.append({
          "WPI_id":wpi,"portwatch_id":portid,
          "original_island_group":p["original_WPI_island_id"],
          "name":p["original_WPI_port_name"],
          "country_source_iso3":p["original_WPI_expected_provider_ISO3"],
          "original_WPI_lat":p["original_WPI_lat"],
          "original_WPI_lon":p["original_WPI_lon"],
          "source_registry_geodesic_km":p["verified_provider_geodesic_km"],
        })
    if len({p["original_island_group"] for p in selected})!=4:
        raise ValueError("Four verified source ports must be in distinct original island groups")
    return selected

def get_source_json(params):
    url=SOURCE_LAYER+"/query?"+urlencode(params)
    if not url.startswith("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/"):
        raise ValueError("Unapproved provider host for WPI4 source sample")
    request=Request(url,headers={"User-Agent":"chocho-20261107-four-WPI-portwatch-edges-bounded-source-probe/1.0",
                                "Accept":"application/json"})
    try:
        with urlopen(request,timeout=30) as r:
            b=r.read(750000+1)
            if len(b)>750000:
                raise ValueError("PortWatch bounded source reply too large")
        obj=json.loads(b.decode("utf-8"))
        if not isinstance(obj,dict) or obj.get("error"):
            return None,{"kind":"ARCGIS_BODY_ERROR","message":str(obj.get("error"))[:280]}
        return obj,None
    except (HTTPError,URLError,TimeoutError,ValueError,ConnectionError) as exc:
        return None,{"kind":type(exc).__name__,"message":str(exc)[:200]}

def run(out_dir):
    source=fixed_samples()
    assigned={x["portwatch_id"]:x for x in source}
    allowed=",".join("'"+x+"'" for x in PROVIDER_IDS)
    predicate="from_id IN ("+allowed+") AND to_id IN ("+allowed+")"
    count,count_error=get_source_json({"where":predicate,"returnCountOnly":"true","f":"json"})
    n=count.get("count") if isinstance(count,dict) else None
    response,row_error=get_source_json({
       "where":predicate,"outFields":",".join(REQUEST_FIELDS),
       "resultRecordCount":str(MAX_SAMPLE_ROWS),
       "returnGeometry":"false","f":"json"
    })
    examples=[x.get("attributes") or {} for x in response.get("features",[])] if isinstance(response,dict) else []
    complete=isinstance(n,int) and n==len(examples) and n<=MAX_SAMPLE_ROWS and count_error is None and row_error is None
    linkrows=[]
    edges={}
    for row in examples:
        a=str(row.get("from_id"));b=str(row.get("to_id"))
        if a not in assigned or b not in assigned:
            raise ValueError("Unexpected PortWatch source port outside four frozen IDs")
        if (str(row.get("from_iso3"))!=assigned[a]["country_source_iso3"] or
            str(row.get("to_iso3"))!=assigned[b]["country_source_iso3"]):
            raise ValueError("PortWatch source edge endpoint country/ISO3 mismatches verified registry source")
        # Sample source has coordinates attached to edges, while provider registry
        # has much more precise canonical lat/lon positions. Check country plus
        # source coordinate proximity to original WPI source, not exact equality.
        for prefix,portid in (("from",a),("to",b)):
            y=assigned[portid]
            phi1=math.radians(float(row[prefix+"_lat"]));phi2=math.radians(float(y["original_WPI_lat"]))
            dl=phi2-phi1;dx=math.radians(float(y["original_WPI_lon"])-float(row[prefix+"_lon"]))
            u=math.sin(dl/2)**2+math.cos(phi1)*math.cos(phi2)*math.sin(dx/2)**2
            km=2*6371.0088*math.asin(math.sqrt(min(1.,max(0.,u))))
            if km>25:
                raise ValueError("Source maritime model edge latitude/longitude not aligned to verified WPI port")
        key=(a,b)
        if key in edges:
            raise ValueError("Provider reported duplicate directed model links, cannot assume binary four-port relation")
        edges[key]=True
        linkrows.append({
          "from_WPI_id":assigned[a]["WPI_id"],
          "to_WPI_id":assigned[b]["WPI_id"],
          "from_WPI_original_island_id":assigned[a]["original_island_group"],
          "to_WPI_original_island_id":assigned[b]["original_island_group"],
          "from_portwatch_id":a,"to_portwatch_id":b,
          "from_name_source":str(row.get("from_name") or "")[:180],
          "to_name_source":str(row.get("to_name") or "")[:180],
          "from_iso3":row.get("from_iso3"),"to_iso3":row.get("to_iso3"),
          "source_maritime_network_distance_field":row.get("distance"),
          "source_maritime_network_transit_days_field":row.get("transit_days"),
          "not_an_individually_observed_vessel_voyage":True,
        })
    pairs=[]
    for a in source:
        for b in source:
            if a["portwatch_id"]==b["portwatch_id"]:
                continue
            pairs.append({
               "from_WPI_id":a["WPI_id"],"to_WPI_id":b["WPI_id"],
               "from_island_id":a["original_island_group"],
               "to_island_id":b["original_island_group"],
               "from_provider_portid":a["portwatch_id"],
               "to_provider_portid":b["portwatch_id"],
               "provider_induced_subgraph_sample_complete":complete,
               "source_has_a_modeled_connection_when_complete":(
                    (a["portwatch_id"],b["portwatch_id"]) in edges if complete else None),
               "provider_record_count_is_not_a_vessel_trip_frequency":True,
            })
    if len(pairs)!=12:
        raise ValueError("Expected exactly 12 possible directed cross-island links in four-port source example")
    unordered=set(tuple(sorted((a,b))) for a,b in edges if a!=b)
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/k).resolve()) for k in ("results","data","provenance")):
        raise ValueError("Never overwrite original WPI source or prior results")
    out.mkdir(parents=True,exist_ok=True)
    receipt={
      "schema":"chocho_WPI4_original_ports_IMF_PortWatch_source_induced_micrograph_v1",
      "status":"FOUR_VERIFIED_WPI_PROVIDER_PORTS_BOUNDED_NETWORK_SOURCE_SAMPLE_COMPLETE",
      "source_layer_URL":SOURCE_LAYER,
      "source_checked_utc":datetime.now(timezone.utc).isoformat(),
      "source_original_publisher":"International Monetary Fund IMF PortWatch, Spillover Simulator",
      "source_provider_full_network_metadata_count_NOT_DOWNLOADED":155521,
      "source_WPI_sample_selection":["28","428","817","1416"],
      "source_verified_WPI_to_PortWatch_nodes":source,
      "source_where_predicate":predicate,
      "provider_4port_induced_subgraph_countOnly":n,
      "provider_sample_feature_rows":len(examples),
      "provider_row_query_cap":MAX_SAMPLE_ROWS,
      "provider_network_query_error":row_error,
      "provider_count_query_error":count_error,
      "source_induced_subgraph_response_complete":complete,
      "source_model_directed_WPI_port_edges":linkrows,
      "source_model_directed_edge_count":len(edges),
      "source_model_unordered_portpair_count":len(unordered),
      "source_all12_cross_original_island_directed_candidate_pairs":pairs,
      "no_full_155521_source_feature_download_performed":True,
      "source_sample_model_edges_NOT_actual_individual_AIS_2020_voyages":True,
      "no_actual_butterfly_or_host_plant_disperal_model":True,
      "source_original_194_WPI_ports_retained_unmodified":True,
    }
    (out/"WPI4_bounded_portwatch_maritime_micrograph_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PORTWATCH_WPI4_INDUCED_GRAPH_SOURCE_JSON="+json.dumps({
       "original_WPI_ports":[x["WPI_id"] for x in source],
       "provider_port_ids":list(PROVIDER_IDS),
       "public_countOnly":n,"sample_returned_rows":len(examples),
       "fully_enumerated_4port_induced_sample":complete,
       "directed_model_links_in_4port_induced_sample":len(edges),
       "unordered_4port_model_link_types":len(unordered),
       "models_not_individual_voyages":True,
       "pair_ids":[[x["from_WPI_id"],x["to_WPI_id"]] for x in linkrows]
    },sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)

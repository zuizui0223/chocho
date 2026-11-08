#!/usr/bin/env python3
"""Audit original source GLSN2018 vs Zenodo2020 OD for PortWatch-confirmed WPI ports.

Port identity: exact union of original 15-of-17 and new 9-of-13 WPI source ports with independent
PortWatch exact-name ISO3 <=25km match and LOCODE. NEVER equate GLSN model
co-service links to real voyages or zeros in one data source to ship absence.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SUB17=ROOT/"results/island51_PortWatch17_geoverification_20261108"
SUB13=ROOT/"results/island51_PortWatch13_geoverification_20261108"
JOIN=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv"
MODEL=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_originalWPI_source_candidate_model_edges.csv"
AIS=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"

def read(p):
    with Path(p).open("rt",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def source_24():
    b=WPI.read_bytes()
    if hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()!=WPI_BLOB:
        raise ValueError("Original WPI194 geographical source Git blob changed")
    wpi=read(WPI)
    port={str(int(r["port_id"])):r for r in wpi}
    if len(wpi)!=194 or len(port)!=194 or len({r["island_id"] for r in wpi})!=51:
        raise ValueError("Original WPI 51 island/194 source port identity altered")
    verified={}
    all_status=[]
    sources=(
        (SUB17,"17",17,15,"ALL17_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED"),
        (SUB13,"13",13,9,"ALL13_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED"),
    )
    for folder,label,total,match_count,status in sources:
        rows=read(folder/("WPI"+label+"_independent_PortWatch_source_geography_LOCODE.csv"))
        independent=json.loads((folder/("independent_WPI"+label+"_source_audit.json")).read_text())
        receipt=json.loads((folder/("WPI"+label+"_PortWatch_source_identity_receipt.json")).read_text())
        if (len(rows)!=total or len({r["original_WPI_port_id"] for r in rows})!=total or
            independent["status"]!=status or receipt["source_verified_WPI_port_count"]!=match_count):
            raise ValueError("Original independent PortWatch WPI"+label+" source identity audit not preserved")
        matched=0
        for r in rows:
            pid=str(int(r["original_WPI_port_id"]))
            if (pid in all_status or pid not in port or
                r["original_island_id"]!=port[pid]["island_id"]):
                raise ValueError("Original source PortWatch identity missing, duplicated or wrong landmass")
            all_status.append(pid)
            if r["original_WPI_independent_geographic_identity_verified"]=="True":
                if (not r["verified_PortWatch_portid"] or not r["verified_PortWatch_LOCODE"] or
                    r["external_registry_box_status"]!="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH" or
                    float(r["original_WPI_to_PortWatch_distance_km"])>25.0):
                    raise ValueError("Unverified source port promoted despite absent exact/ISO/geographic status")
                verified[pid]=r
                matched+=1
            else:
                if r["verified_PortWatch_portid"] or r["verified_PortWatch_LOCODE"]:
                    raise ValueError("Unverified original WPI candidate has a provider ID or LOCODE")
        if matched!=match_count:
            raise ValueError("Source WPI"+label+" full cohort verified-match count differs")
    if len(all_status)!=30 or len(verified)!=24 or len({r["original_island_id"] for r in verified.values()})!=17:
        raise ValueError("Merged source population must be 30 candidates, 24 verified ports and 17 islands")
    if len({r["verified_PortWatch_portid"] for r in verified.values()})!=24:
        raise ValueError("One external PortWatch provider port ID assigned to different original WPI ports")
    return verified,port

def to_csv(p,records,fields):
    with Path(p).open("wt",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

def run(out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/s).resolve()) for s in ("data","results","provenance")):
        raise ValueError("Do not overwrite original scientific results")
    out.mkdir(parents=True,exist_ok=True)
    verified,wpi=source_24()
    history=read(JOIN)
    hist_byid={str(int(r["original_WPI_port_id"])):r for r in history}
    glsn=read(MODEL)
    ais=read(AIS)
    if len(history)!=194 or len(hist_byid)!=194 or len(glsn)!=37 or len(ais)!=94:
        raise ValueError("Original complete GLSN/AIS/UNCTAD source cohorts have changed")
    model_edges={}
    for edge in glsn:
        a,b=(str(int(edge[k])) for k in ("original_WPI_candidate_port_1","original_WPI_candidate_port_2"))
        key=tuple(sorted((a,b),key=int))
        if key in model_edges:
            raise ValueError("GLSN original undirected candidate model edge repeated")
        model_edges[key]=edge
    directed_ais={}
    for edge in ais:
        a,b=(str(int(edge[k])) for k in ("origin_WPI_port_id","dest_WPI_port_id"))
        if a==b or (a,b) in directed_ais:
            raise ValueError("AIS source port OD has repeated directed original port link")
        n=int(edge["source_observed_2020_triplet_count"])
        if n<=0:
            raise ValueError("Zero/invalid original AIS2020 source triplet in reported edge")
        directed_ais[(a,b)]=n
    ids=sorted(verified,key=int)
    pairs=[]
    model_source_island=Counter()
    ais_source_island=Counter()
    ais_direct_island=Counter()
    for a,b in itertools.combinations(ids,2):
        x,y=verified[a],verified[b]
        ia,ib=x["original_island_id"],y["original_island_id"]
        key=(a,b)
        model=model_edges.get(key)
        abin=directed_ais.get((a,b),0)
        bain=directed_ais.get((b,a),0)
        both=abin+bain
        model_present=model is not None
        if ia!=ib:
            ip=tuple(sorted((ia,ib),key=int))
            if model_present:model_source_island[ip]+=1
            if both:ais_source_island[ip]+=both
            if abin:ais_direct_island[(ia,ib)]+=abin
            if bain:ais_direct_island[(ib,ia)]+=bain
        ra=hist_byid[a];rb=hist_byid[b]
        for pid,r in ((a,ra),(b,rb)):
            if (r["GLSN2018_status"]!="UNIQUE_ISO3_TO_ISO2_AND_NAME_CANDIDATE_NOT_GEO_VERIFIED" or
                not r["UNCTAD_PLSCIQ_match_status"].startswith("UNIQUE_")):
                raise ValueError("Independent geography verified source port not in original GLSN+UNCTAD double-source panel: "+pid)
        def history_value(r):
            return (r["UNCTAD_2006_2010_port_index_mean"] if r["UNCTAD_2006_2010_available_quarters"] and
                    int(r["UNCTAD_2006_2010_available_quarters"])>0 else "")
        pairs.append({
            "original_WPI_port_a":a,"original_WPI_port_b":b,
            "original_island_a":ia,"original_island_b":ib,
            "portwatch_verified_portid_a":x["verified_PortWatch_portid"],
            "portwatch_verified_portid_b":y["verified_PortWatch_portid"],
            "portwatch_LOCODE_a":x["verified_PortWatch_LOCODE"],
            "portwatch_LOCODE_b":y["verified_PortWatch_LOCODE"],
            "both_original_ports_independently_geo_verified":True,
            "GLSN2018_model_co_service_edge":model_present,
            "GLSN2018_source_GC_km":model["original_source_GC_km"] if model else "",
            "GLSN2018_source_weight_MODEL_NOT_SHIP_TRIPS":model["original_source_weight_NOT_TRAFFIC_VOLUME"] if model else "",
            "Zenodo2020_AIS_OD_source_rows_a_to_b":abin,
            "Zenodo2020_AIS_OD_source_rows_b_to_a":bain,
            "Zenodo2020_AIS_OD_source_rows_both":both,
            "both_sources_record_port_pair":model_present and both>0,
            "UNCTAD_2006_2010_available_quarters_a":ra["UNCTAD_2006_2010_available_quarters"],
            "UNCTAD_2006_2010_available_quarters_b":rb["UNCTAD_2006_2010_available_quarters"],
            "UNCTAD_2006_2010_index_mean_a":history_value(ra),
            "UNCTAD_2006_2010_index_mean_b":history_value(rb),
            "GLSN_model_edge_absence_is_NOT_no_shipping":True,
            "AIS_source_zero_is_NOT_actual_no_shipping":True,
            "no_causal_butterfly_or_shipborne_host_inference":True,
        })
    if len(pairs)!=276:
        raise ValueError("Original 24 port choose 2 original WPI verified pairs missing")
    source_port_count=Counter(x["original_island_id"] for x in verified.values())
    island_ids=sorted({w["island_id"] for w in wpi.values()},key=int)
    source_island_pairs=[]
    for i,ia in enumerate(island_ids):
        for ib in island_ids[i+1:]:
            pair=(ia,ib)
            source_island_pairs.append({
                "original_island_a":ia,"original_island_b":ib,
                "independently_geo_verified_source_ports_on_a":source_port_count[ia],
                "independently_geo_verified_source_ports_on_b":source_port_count[ib],
                "GLSN2018_model_co_service_edges_between_verified_ports":model_source_island[pair],
                "Zenodo2020_source_trips_between_geo_verified_ports":ais_source_island[pair],
                "GLSN_and_Zenodo_both_have_source_pair_type":model_source_island[pair]>0 and ais_source_island[pair]>0,
                "unverified_port_or_absent_source_link_NOT_ship_traffic_absence":True,
                "actual_shipborne_butterfly_or_host_cargo_UNVERIFIED":True,
            })
    if len(source_island_pairs)!=1275:
        raise ValueError("Original 51-island unordered pair ledger incomplete")
    to_csv(out/"WPI24_geoverified_all276_GLSN_UNCTAD_AIS_port_pair_comparison.csv",pairs,list(pairs[0]))
    to_csv(out/"original51_all1275_geo_verified24_source_island_pair_observation.csv",
           source_island_pairs,list(source_island_pairs[0]))
    common=[r for r in pairs if r["both_sources_record_port_pair"]]
    hist16=sum(int(hist_byid[pid]["UNCTAD_2006_2010_available_quarters"] or 0)>=16 for pid in ids)
    rec={
       "schema":"chocho_source_only_IMF_geoverified24_port_GLSN_UNCTAD_AIS_network_support_v1",
       "status":"ALL_276_GEO_VERIFIED_PORT_PAIRS_AND_ORIGINAL_1275_ISLAND_PAIRS_RECONSTRUCTED",
       "original_WPI_ports":194,"original_islands":51,
       "source_selected_WPI_ports":30,
       "independently_geo_verified_WPI_ports":len(ids),
       "independently_geo_verified_original_islands":len(source_port_count),
       "all_24_source_port_pairs":len(pairs),
       "model_GLSN2018_co_service_port_pairs":sum(bool(r["GLSN2018_model_co_service_edge"]) for r in pairs),
       "AIS2020_original_source_undirected_port_pairs":sum(r["Zenodo2020_AIS_OD_source_rows_both"]>0 for r in pairs),
       "AIS2020_original_source_directed_port_links":sum(int(r["Zenodo2020_AIS_OD_source_rows_a_to_b"]>0)+int(r["Zenodo2020_AIS_OD_source_rows_b_to_a"]>0) for r in pairs),
       "AIS2020_source_inferred_OD_episode_rows":sum(r["Zenodo2020_AIS_OD_source_rows_both"] for r in pairs),
       "modelled_and_AIS_source_shared_port_pair_count":len(common),
       "modelled_and_AIS_source_shared_original_WPI_port_pairs":[[x["original_WPI_port_a"],x["original_WPI_port_b"]] for x in common],
       "GLSN2018_original_island_pair_types":len(model_source_island),
       "AIS2020_original_island_pair_types":len(ais_source_island),
       "source_both_GLSN_and_AIS_original_island_pair_types":sum(ia in ais_source_island for ia in model_source_island),
       "UNCTAD2006_2010_16of20_source_ports":hist16,
       "UNCTAD2006_2010_16of20_missing_index_ports":len(ids)-hist16,
       "WPI24_source_WPI_identities":[{"wpi_id":x,"original_island":verified[x]["original_island_id"],
            "portwatch_id":verified[x]["verified_PortWatch_portid"],
            "locode":verified[x]["verified_PortWatch_LOCODE"]} for x in ids],
       "full_194_WPI17_port_identity_receipt_sha256":sha(SUB17/"WPI17_PortWatch_source_identity_receipt.json"),
       "full_194_WPI13_port_identity_receipt_sha256":sha(SUB13/"WPI13_PortWatch_source_identity_receipt.json"),
       "original_source_GLSN_pair_file_SHA256":sha(MODEL),
       "original_source_2020_OD_file_SHA256":sha(AIS),
       "original_historical_UNCTAD_194_file_SHA256":sha(JOIN),
       "outputs_sha256":{
           "port_pairs":sha(out/"WPI24_geoverified_all276_GLSN_UNCTAD_AIS_port_pair_comparison.csv"),
           "island_pairs":sha(out/"original51_all1275_geo_verified24_source_island_pair_observation.csv"),
       },
       "source_only_eligible_cohort_not_representative_all_original_ports":True,
       "GLSN_source_edges_model_co_membership_NOT_observed_ship_OD":True,
       "UNCTAD_quarterly_index_NOT_real_shipping_volumes_or_OD":True,
       "Zenodo_source_OD_2020_not_2001_2010_shipping":True,
       "no_imputed_shipping_absence_or_causal_butterfly_claim":True,
       "computed_utc":datetime.now(timezone.utc).isoformat(),
    }
    (out/"WPI24_3_source_verified_port_graph_audit.json").write_text(
        json.dumps(rec,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PORTWATCH24_CROSS_SOURCE_GRAPH_RESULT_JSON="+json.dumps({
        k:rec[k] for k in ("status","independently_geo_verified_WPI_ports",
            "independently_geo_verified_original_islands","all_24_source_port_pairs",
            "model_GLSN2018_co_service_port_pairs",
            "AIS2020_original_source_undirected_port_pairs",
            "AIS2020_original_source_directed_port_links",
            "AIS2020_source_inferred_OD_episode_rows",
            "modelled_and_AIS_source_shared_port_pair_count",
            "modelled_and_AIS_source_shared_original_WPI_port_pairs",
            "GLSN2018_original_island_pair_types","AIS2020_original_island_pair_types",
            "source_both_GLSN_and_AIS_original_island_pair_types","UNCTAD2006_2010_16of20_source_ports")},sort_keys=True),flush=True)
    return rec

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    run(p.parse_args().out_dir)

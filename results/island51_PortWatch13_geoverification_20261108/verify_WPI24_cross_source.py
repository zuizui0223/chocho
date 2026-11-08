#!/usr/bin/env python3
"""Independent source-data audit of all 276 verified-port and 1,275 original island pairs.

No fitted ecological responses. Reconstruct exact 2018 source MODEL co-service
edges, 2006-2010 indexed UNCTAD port-time source, and separately recorded 2020
AIS-inferred port OD episodes, using only original 24 geoverified WPI port IDs.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,itertools,math
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
P17=ROOT/"results/island51_PortWatch17_geoverification_20261108"
P13=ROOT/"results/island51_PortWatch13_geoverification_20261108"
GL=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_originalWPI_source_candidate_model_edges.csv"
JOIN=ROOT/"results/island51_maritime_raw_recovery_20261108/GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv"
AIS=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv"
ORIG=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
NAME="WPI24_geoverified_all276_GLSN_UNCTAD_AIS_port_pair_comparison.csv"
ISLAND="original51_all1275_geo_verified24_source_island_pair_observation.csv"

def read(path):
    with Path(path).open("rt",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def key(a,b):
    return tuple(sorted((str(int(a)),str(int(b))),key=int))

def compare_float(a,b,label):
    if not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError("Changed original source "+label)

def audit(directory):
    path=Path(directory)
    contents=ORIG.read_bytes()
    if hashlib.sha1(b"blob "+str(len(contents)).encode()+b"\0"+contents).hexdigest()!=BLOB:
        raise ValueError("Original WPI source changed")
    original=read(ORIG)
    if len(original)!=194 or len({r["island_id"] for r in original})!=51:
        raise ValueError("Original 51 islands / 194 WPI ports changed")
    wpi={str(int(x["port_id"])):x for x in original}
    old=read(P17/"WPI17_independent_PortWatch_source_geography_LOCODE.csv")
    new=read(P13/"WPI13_independent_PortWatch_source_geography_LOCODE.csv")
    r17=json.loads((P17/"independent_WPI17_source_audit.json").read_text())
    r13=json.loads((P13/"independent_WPI13_source_audit.json").read_text())
    if (len(old)!=17 or len(new)!=13 or
        r17["status"]!="ALL17_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED" or
        r13["status"]!="ALL13_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED"):
        raise ValueError("Source-locked old/new independent PortWatch registry audits missing")
    all_candidates=old+new
    if len({x["original_WPI_port_id"] for x in all_candidates})!=30:
        raise ValueError("Expected exactly 30 original two-provider WPI source candidates")
    ids={str(int(x["original_WPI_port_id"])):x for x in all_candidates if x["verified_PortWatch_portid"]}
    if len(ids)!=24 or len({x["original_island_id"] for x in ids.values()})!=17:
        raise ValueError("Independent PortWatch source exactly 24 original ports in 17 islands required")
    if len({x["verified_PortWatch_portid"] for x in ids.values()})!=24:
        raise ValueError("Multiple WPI ports share one PortWatch provider ID")
    for r in all_candidates:
        pid=str(int(r["original_WPI_port_id"]))
        if pid not in wpi or r["original_island_id"]!=wpi[pid]["island_id"]:
            raise ValueError("Original WPI port/island source identity changed")
        if pid not in ids and (r["verified_PortWatch_portid"] or
                              r["original_WPI_independent_geographic_identity_verified"].lower()!="false"):
            raise ValueError("Unverified prior/new provider port incorrectly promoted")
    historical=read(JOIN)
    history={str(int(x["original_WPI_port_id"])):x for x in historical}
    glsn=read(GL)
    ais=read(AIS)
    if len(glsn)!=37 or len(ais)!=94 or len(history)!=194:
        raise ValueError("Original network/UNCTAD/AIS source panel changed")
    model={}
    for edge in glsn:
        p=key(edge["original_WPI_candidate_port_1"],edge["original_WPI_candidate_port_2"])
        if p in model:
            raise ValueError("Duplicate modelled 2018 source port pair")
        model[p]=edge
    observed={}
    for row in ais:
        a,b=(str(int(row[k])) for k in ("origin_WPI_port_id","dest_WPI_port_id"))
        if a==b or (a,b) in observed:
            raise ValueError("Original 2020 AIS source directed port edge invalid/duplicated")
        observed[(a,b)]=int(row["source_observed_2020_triplet_count"])
        if observed[(a,b)]<=0:raise ValueError("Original source AIS episode count nonpositive")
    table=read(path/NAME)
    if len(table)!=276:
        raise ValueError("Require all 276 unordered verified original port pair combinations")
    seen=set()
    island_model=Counter();island_ais=Counter()
    model_count=0;ais_count=0;directed_count=0;ais_episode_rows=0;overlap=[]
    for row in table:
        a,b=row["original_WPI_port_a"],row["original_WPI_port_b"]
        pair=key(a,b)
        if pair in seen or pair!=(a,b) or a not in ids or b not in ids:
            raise ValueError("Invalid or repeated source WPI24 original port pair")
        seen.add(pair)
        ia,ib=ids[a]["original_island_id"],ids[b]["original_island_id"]
        if (row["original_island_a"]!=ia or row["original_island_b"]!=ib or
            row["portwatch_verified_portid_a"]!=ids[a]["verified_PortWatch_portid"] or
            row["portwatch_verified_portid_b"]!=ids[b]["verified_PortWatch_portid"] or
            row["portwatch_LOCODE_a"]!=ids[a]["verified_PortWatch_LOCODE"] or
            row["portwatch_LOCODE_b"]!=ids[b]["verified_PortWatch_LOCODE"]):
            raise ValueError("Published original WPI/PortWatch ID/LOCODE identity differs")
        gl= model.get(pair)
        source_modelled=bool(gl)
        if (row["GLSN2018_model_co_service_edge"].lower()=="true")!=source_modelled:
            raise ValueError("Original 2018 source co-service edge identity changed")
        if gl:
            compare_float(row["GLSN2018_source_GC_km"],gl["original_source_GC_km"],"GLSN source GC")
            compare_float(row["GLSN2018_source_weight_MODEL_NOT_SHIP_TRIPS"],gl["original_source_weight_NOT_TRAFFIC_VOLUME"],"GLSN model weight")
            model_count+=1
        elif row["GLSN2018_source_GC_km"] or row["GLSN2018_source_weight_MODEL_NOT_SHIP_TRIPS"]:
            raise ValueError("A nonexistent source model edge was assigned made-up GC/weight")
        ab=observed.get((a,b),0);ba=observed.get((b,a),0)
        if int(row["Zenodo2020_AIS_OD_source_rows_a_to_b"])!=ab or int(row["Zenodo2020_AIS_OD_source_rows_b_to_a"])!=ba or int(row["Zenodo2020_AIS_OD_source_rows_both"])!=ab+ba:
            raise ValueError("Original inferred 2020 source AIS port direction/episode counts changed")
        ais_count+=int(ab+ba>0)
        directed_count+=int(ab>0)+int(ba>0)
        ais_episode_rows+=ab+ba
        if (row["both_sources_record_port_pair"].lower()=="true")!=(source_modelled and ab+ba>0):
            raise ValueError("Original source cross-network port pair intersection changed")
        if source_modelled and ab+ba>0:overlap.append(pair)
        if ia!=ib:
            ip=key(ia,ib)
            if source_modelled:island_model[ip]+=1
            island_ais[ip]+=ab+ba
        ha,hb=history[a],history[b]
        for side,original_port in (("a",ha),("b",hb)):
            field="UNCTAD_2006_2010_available_quarters_"+side
            if row[field]!=original_port["UNCTAD_2006_2010_available_quarters"]:
                raise ValueError("Original historical port index quarter support changed")
            source_mean=original_port["UNCTAD_2006_2010_port_index_mean"]
            saved_mean=row["UNCTAD_2006_2010_index_mean_"+side]
            if source_mean:
                compare_float(saved_mean,source_mean,"UNCTAD historical index source")
            elif saved_mean:
                raise ValueError("Missing historical UNCTAD source index was converted to zero")
        for flag in ("GLSN_model_edge_absence_is_NOT_no_shipping",
                     "AIS_source_zero_is_NOT_actual_no_shipping",
                     "no_causal_butterfly_or_shipborne_host_inference"):
            if row[flag].lower()!="true":
                raise ValueError("False source-ship absence or causal label: "+flag)
    if len(seen)!=276 or seen!=set(itertools.combinations(sorted(ids,key=int),2)):
        raise ValueError("Not every original 24-port pair audited exactly once")
    isTable=read(path/ISLAND)
    if len(isTable)!=1275:
        raise ValueError("Original all 1275 island-pair source comparison population omitted")
    cnt=Counter(r["original_island_id"] for r in ids.values())
    all_island=sorted({r["island_id"] for r in original},key=int)
    expected_island_pairs=set(itertools.combinations(all_island,2))
    checked=set()
    for row in isTable:
        pair=(row["original_island_a"],row["original_island_b"])
        if pair not in expected_island_pairs or pair in checked:
            raise ValueError("Duplicated/nonoriginal 51-island pair source")
        checked.add(pair)
        if (int(row["independently_geo_verified_source_ports_on_a"])!=cnt[pair[0]] or
            int(row["independently_geo_verified_source_ports_on_b"])!=cnt[pair[1]] or
            int(row["GLSN2018_model_co_service_edges_between_verified_ports"])!=island_model[pair] or
            int(row["Zenodo2020_source_trips_between_geo_verified_ports"])!=island_ais[pair] or
            (row["GLSN_and_Zenodo_both_have_source_pair_type"].lower()=="true")!=(
               island_model[pair]>0 and island_ais[pair]>0)):
            raise ValueError("Original source island-pair comparisons differ")
        if row["unverified_port_or_absent_source_link_NOT_ship_traffic_absence"].lower()!="true":
            raise ValueError("Source observation absence illegally recoded as biological no ship")
    if len(checked)!=1275:
        raise ValueError("Original source all 1275 island pairs incomplete")
    # The original PR147 15-port panel is a strict immutable subset of new 24-port outcomes.
    prior_rows=read(P17/"WPI15_geoverified_all105_GLSN_UNCTAD_AIS_port_pair_comparison.csv")
    original_fifteen={str(int(p["original_WPI_port_id"])) for p in old if p["verified_PortWatch_portid"]}
    if len(original_fifteen)!=15 or len(prior_rows)!=105:
        raise ValueError("Original PR147 15-port/105-pair source outcomes unavailable")
    parent_pairs={key(x["original_WPI_port_a"],x["original_WPI_port_b"]):x for x in prior_rows}
    child_pairs={key(x["original_WPI_port_a"],x["original_WPI_port_b"]):x
                 for x in table if x["original_WPI_port_a"] in original_fifteen
                                 and x["original_WPI_port_b"] in original_fifteen}
    if len(parent_pairs)!=105 or len(child_pairs)!=105 or set(parent_pairs)!=set(child_pairs):
        raise ValueError("Original 15-port study not preserved as exact subset")
    for pair in parent_pairs:
        previous,current=parent_pairs[pair],child_pairs[pair]
        for col in ("GLSN2018_model_co_service_edge","GLSN2018_source_GC_km",
                    "GLSN2018_source_weight_MODEL_NOT_SHIP_TRIPS",
                    "Zenodo2020_AIS_OD_source_rows_a_to_b",
                    "Zenodo2020_AIS_OD_source_rows_b_to_a",
                    "Zenodo2020_AIS_OD_source_rows_both",
                    "UNCTAD_2006_2010_available_quarters_a",
                    "UNCTAD_2006_2010_available_quarters_b",
                    "UNCTAD_2006_2010_index_mean_a",
                    "UNCTAD_2006_2010_index_mean_b"):
            if previous[col]!=current[col]:
                raise ValueError("Previously independently audited 15-port source value changed: "+col)
    receipt=json.loads((path/"WPI24_3_source_verified_port_graph_audit.json").read_text())
    derived={
      "independently_geo_verified_WPI_ports":len(ids),
      "independently_geo_verified_original_islands":len(cnt),
      "all_15_source_port_pairs":len(seen),
      "model_GLSN2018_co_service_port_pairs":model_count,
      "AIS2020_original_source_undirected_port_pairs":ais_count,
      "AIS2020_original_source_directed_port_links":directed_count,
      "AIS2020_source_inferred_OD_episode_rows":ais_episode_rows,
      "modelled_and_AIS_source_shared_port_pair_count":len(overlap),
      "GLSN2018_original_island_pair_types":len(island_model),
      "AIS2020_original_island_pair_types":sum(x>0 for x in island_ais.values()),
      "source_both_GLSN_and_AIS_original_island_pair_types":sum(x>0 and island_ais[k]>0 for k,x in island_model.items()),
      "UNCTAD2006_2010_16of20_source_ports":sum(int(history[i]["UNCTAD_2006_2010_available_quarters"] or 0)>=16 for i in ids)
    }
    for k,v in derived.items():
        if receipt[k]!=v:
            raise ValueError("Published source-only graph/port result receipt does not reproduce: "+k)
    if receipt["modelled_and_AIS_source_shared_original_WPI_port_pairs"]!=[list(x) for x in overlap]:
        raise ValueError("Original GLSN versus 2020 OD pair identities were replaced")
    result={"schema":"chocho_independently_geo_verified_WPI24_3marine_source_audit_v1",
       "status":"ALL_276_VERIFIED_PORT_PAIRS_AND_1275_ORIGINAL_ISLAND_PAIRS_RECONSTRUCTED",
       **derived,
       "source_modelled_GLSN_not_actual_2020_OD":True,
       "index_missing_not_zero":True,
       "no_194_port_source_island_traffic_absences_inferred":True}
    (path/"independent_WPI24_3source_graph_audit.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("PORTWATCH15_CROSS_SOURCE_INDEPENDENT_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--results-dir",required=True,type=Path)
    audit(p.parse_args().results_dir)

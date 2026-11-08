#!/usr/bin/env python3
"""Audit complete GLSN2018 source graph and candidate match to original 194 WPI ports.

GLSN study source: MIT-licensed original Excel, archived as full extracted sheets
after publisher Git blob verification. The 2018 edge list is liner-service graph
co-membership, not a verified consecutive physical vessel journey.
Original GLSN port sheet HAS ISO3 but NO coordinates: country/name joins are
candidate identities only, not verified geolocation/port-object mappings.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import pycountry

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/island51_maritime_raw_recovery_20261108"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
GLSNPORT=SOURCE/"GLSN2018_ports_original_sheet_1.csv"
GLSNEDGE=SOURCE/"GLSN2018_edges_original_sheet_1.csv"
GLSNRECEIPT=SOURCE/"GLSN2018_recovered_source_receipt.json"
UNCTAD=SOURCE/"UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv"
AIS=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/original194_source_matched_port_coverage.csv"

def read(path):
    with Path(path).open("rt",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def write(path,rows,fields):
    with Path(path).open("wt",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def gitblob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def canonical(s):
    s=unicodedata.normalize("NFKD",str(s)).encode("ascii","ignore").decode("ascii").upper()
    tokens=re.findall(r"[A-Z0-9]+",s)
    if not tokens:return ""
    if len(tokens)>2 and tokens[:2]==["PORT","OF"]:
        tokens=tokens[2:]
    elif len(tokens)>1 and tokens[0]=="PORT":
        tokens=tokens[1:]
    while len(tokens)>1 and tokens[-1] in {"PORT","HARBOR","HARBOUR","TERMINAL"}:
        tokens.pop()
    if not tokens or tokens in (["PORT"],["NAN"]):
        return ""
    return "".join(tokens)

def iso3_to_iso2(value):
    result=pycountry.countries.get(alpha_3=str(value).strip().upper())
    return result.alpha_2 if result else ""

def verify_source():
    if gitblob(WPI)!=WPI_BLOB:
        raise ValueError("Original WPI194 source island/port geography changed")
    glsn_meta=json.loads(GLSNRECEIPT.read_text())
    if glsn_meta.get("status")!="ORIGINAL_PUBLISHER_GIT_BLOB_XLSX_AND_FULL_EXTRACTED_TABLES_VERIFIED":
        raise ValueError("Publisher GLSN original XLSX source not acquired/verified")
    source_specs={x["role"]:x for x in glsn_meta["sources"]}
    for role,file in (("PORTS",GLSNPORT),("EDGES",GLSNEDGE)):
        sheets=source_specs[role]["sheets"]
        if len(sheets)!=1 or sheets[0]["extracted_csv_sha256"]!=digest(file):
            raise ValueError("Derived original GLSN spreadsheet data do not match publisher-pinned extraction")
    ports=read(GLSNPORT)
    edges=read(GLSNEDGE)
    original=read(WPI)
    unctad=read(UNCTAD)
    ais=read(AIS)
    if (len(ports)!=887 or len(edges)!=28104 or len(original)!=194 or
        len(unctad)!=194 or len(ais)!=194 or
        len({p["island_id"] for p in original})!=51):
        raise ValueError("The original 887 GLSN ports/28104 edges/194 WPI sources changed")
    if sum(x["included"].lower()=="true" for x in original)!=182:
        raise ValueError("Original 182 source marine qualified ports changed")
    if not set(("port","ISO3","Country Name")).issubset(ports[0]) or not set(("port1","port2","GC","weight")).issubset(edges[0]):
        raise ValueError("Original publisher GLSN port and edge schema changed")
    return ports,edges,original,unctad,ais,glsn_meta

def run(outdir):
    out=Path(outdir).resolve()
    if any(out.is_relative_to((ROOT/s).resolve()) for s in ("data","results","provenance")):
        raise ValueError("Never overwrite original source/legacy scientific results")
    out.mkdir(parents=True,exist_ok=True)
    ports,edges,original,unctad,ais,source_meta=verify_source()
    # Globally same-named ports in GLSN EDGE list are not disambiguated by country,
    # so only a unique normalized provider name can be promoted to a network node.
    provider_name_count=Counter(canonical(p["port"]) for p in ports)
    provider_name_country_count=Counter((iso3_to_iso2(p["ISO3"]),canonical(p["port"])) for p in ports)
    original_name_country_count=Counter((p["country_code"],canonical(p["port_name"])) for p in original)
    provider_by_key=defaultdict(list)
    unresolved_country=Counter()
    for p in ports:
        iso2=iso3_to_iso2(p["ISO3"])
        if not iso2:
            unresolved_country[p["ISO3"]]+=1
        provider_by_key[(iso2,canonical(p["port"]))].append(p)
    unctad_by_id={str(int(p["original_WPI_id"])):p for p in unctad}
    ais_by_id={str(int(p["original_WPI_port_id"])):p for p in ais}
    if len(unctad_by_id)!=194 or len(ais_by_id)!=194:
        raise ValueError("Original UNCTAD or true AIS WPI port identity ledger duplicated")
    matches=[]
    glsn_name_to_wpi={}
    for p in original:
        pid=str(int(p["port_id"]))
        identity=(p["country_code"],canonical(p["port_name"]))
        hits=provider_by_key[identity]
        nameglobally_unique=len(hits)==1 and provider_name_count[canonical(hits[0]["port"])]==1
        strong=(original_name_country_count[identity]==1 and nameglobally_unique)
        hit=hits[0] if strong else None
        if strong:
            provider=hit["port"]
            if provider in glsn_name_to_wpi:
                raise ValueError("A single GLSN source name linked to two distinct WPI source ports")
            glsn_name_to_wpi[provider]=pid
        old=unctad_by_id[pid]
        aisrow=ais_by_id[pid]
        if (str(old["original_island_id"])!=str(p["island_id"]) or
            str(aisrow["island_id"])!=str(p["island_id"])):
            raise ValueError("Matched source original island ID changed across independent panels")
        matches.append({
            "original_WPI_port_id":pid,"island_id":p["island_id"],
            "original_WPI_country_ISO2":p["country_code"],
            "original_WPI_name":p["port_name"],
            "GLSN2018_status":("UNIQUE_ISO3_TO_ISO2_AND_NAME_CANDIDATE_NOT_GEO_VERIFIED" if strong
                else "PROVIDER_NAME_AMBIGUOUS_FOR_EDGE_IDS" if len(hits)>=1 else "UNMATCHED"),
            "GLSN2018_provider_port_name":hit["port"] if hit else "",
            "GLSN2018_provider_ISO3":hit["ISO3"] if hit else "",
            "GLSN2018_provider_country":hit["Country Name"] if hit else "",
            "GLSN2018_original_coords_INDEPENDENTLY_VERIFIED":False,
            "UNCTAD_PLSCIQ_match_status":old["status"],
            "UNCTAD_2006_2010_available_quarters":old["available_2006_2010_quarters"],
            "UNCTAD_2006_2010_port_index_mean":old["mean_2006_2010_PLCSIQ_old_index_base"],
            "Zenodo2020_actual_WPI_port_identity_verified":str(aisrow["verified_external_2020_endpoint"]).lower()=="true",
            "source_data_roles":"GLSN_modelled_liner_network_UNCTAD_port_connectivity_index_Zenodo2020_AIS_inferred_OD",
        })
    # Verify complete weighted network's directional duplication; preserve provenance.
    all_names={p["port"] for p in ports}
    exact_pair_rows=defaultdict(list)
    unknown_endpoints=Counter()
    invalid_weight=0
    bad_GC=0
    self_loops=0
    for r in edges:
        a,b=r["port1"],r["port2"]
        if a not in all_names:unknown_endpoints[a]+=1
        if b not in all_names:unknown_endpoints[b]+=1
        if a==b:self_loops+=1
        try:
            gc=float(r["GC"]);weight=float(r["weight"])
            if not math.isfinite(weight) or weight<0:invalid_weight+=1
            if not math.isfinite(gc) or gc<0:bad_GC+=1
        except (ValueError,TypeError):
            invalid_weight+=1;bad_GC+=1
        exact_pair_rows[tuple(sorted((a,b)))].append(r)
    undirected=len(exact_pair_rows)
    symmetric_name_pairs=sum(len(rs)==2 for rs in exact_pair_rows.values())
    asymmetric_counts=undirected-symmetric_name_pairs
    mismatched_reverse=0
    for pair,rs in exact_pair_rows.items():
        if len(rs)==2:
            a,b=pair
            if {(r["port1"],r["port2"]) for r in rs}!={(a,b),(b,a)}:
                mismatched_reverse+=1
            else:
                weights=[float(r["weight"]) for r in rs]
                gcs=[float(r["GC"]) for r in rs]
                if abs(weights[0]-weights[1])>1e-9 or abs(gcs[0]-gcs[1])>1e-9:
                    mismatched_reverse+=1
    if invalid_weight or bad_GC:
        raise ValueError("GLSN original source edges include invalid weight/GC numerical measurements")
    # Do not pretend source undirected model edges are direct ship voyages.
    candidate_edges=[]
    modelled_between_island=Counter()
    total_modelled_weight=defaultdict(float)
    raw_directed_in_candidate=0
    for pair,rs in sorted(exact_pair_rows.items()):
        a,b=pair
        if a in glsn_name_to_wpi and b in glsn_name_to_wpi and a!=b:
            wpa,wpb=glsn_name_to_wpi[a],glsn_name_to_wpi[b]
            if wpa==wpb:continue
            raw_directed_in_candidate+=len(rs)
            # Original 2018 source may have mirrored rows; one unordered edge only.
            gc=float(rs[0]["GC"]);weight=float(rs[0]["weight"])
            pa=next(p for p in original if str(int(p["port_id"]))==wpa)
            pb=next(p for p in original if str(int(p["port_id"]))==wpb)
            ia,ib=pa["island_id"],pb["island_id"]
            if ia!=ib:
                key=tuple(sorted((ia,ib),key=int))
                modelled_between_island[key]+=1
                total_modelled_weight[key]+=weight
            candidate_edges.append({
                "GLSN_provider_port_1":a,"GLSN_provider_port_2":b,
                "original_WPI_candidate_port_1":wpa,"original_WPI_candidate_port_2":wpb,
                "original_island_1":ia,"original_island_2":ib,
                "original_source_GC_km":gc,"original_source_weight_NOT_TRAFFIC_VOLUME":weight,
                "GLSN_original_directed_rows_for_undirected_pair":len(rs),
                "port_name_country_uniqueness_only_NOT_geographically_verified":True,
                "edge_is_liner_service_co_membership_not_actual_direct_voyage":True,
            })
    # Full original 1275 island-pair universe, missing source edges NOT traffic absence.
    source_wpi_by_island=Counter(r["island_id"] for r in matches
        if r["GLSN2018_status"]=="UNIQUE_ISO3_TO_ISO2_AND_NAME_CANDIDATE_NOT_GEO_VERIFIED")
    pair_rows=[]
    ids=sorted({p["island_id"] for p in original},key=int)
    for i,ia in enumerate(ids):
        for ib in ids[i+1:]:
            key=(ia,ib)
            pair_rows.append({
                "original_island_a":ia,"original_island_b":ib,
                "country_name_candidate_ports_island_a":source_wpi_by_island[ia],
                "country_name_candidate_ports_island_b":source_wpi_by_island[ib],
                "GLSN2018_modelled_undirected_copresence_edges_in_matched_candidate_ports":modelled_between_island[key],
                "GLSN2018_modelled_source_weight_sum_NOT_vessel_trip_count":total_modelled_weight[key],
                "missing_model_source_edge_is_not_absence_of_ship_traffic":True,
                "original_source_island_pair_geographical_identity_NOT_verified":True,
            })
    if len(pair_rows)!=1275:
        raise ValueError("Original source 51 island-pair universe incomplete")
    portfields=list(matches[0])
    write(out/"GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv",matches,portfields)
    if candidate_edges:
        write(out/"GLSN2018_originalWPI_source_candidate_model_edges.csv",candidate_edges,list(candidate_edges[0]))
    else:
        write(out/"GLSN2018_originalWPI_source_candidate_model_edges.csv",[],
              ["GLSN_provider_port_1","GLSN_provider_port_2","original_WPI_candidate_port_1",
               "original_WPI_candidate_port_2","original_island_1","original_island_2",
               "original_source_GC_km","original_source_weight_NOT_TRAFFIC_VOLUME",
               "GLSN_original_directed_rows_for_undirected_pair",
               "port_name_country_uniqueness_only_NOT_geographically_verified",
               "edge_is_liner_service_co_membership_not_actual_direct_voyage"])
    write(out/"GLSN2018_original51_all1275_modelled_island_pair_candidates.csv",pair_rows,list(pair_rows[0]))
    cand=[x for x in matches if x["GLSN2018_status"]=="UNIQUE_ISO3_TO_ISO2_AND_NAME_CANDIDATE_NOT_GEO_VERIFIED"]
    unct=[x for x in matches if x["UNCTAD_PLSCIQ_match_status"]=="UNIQUE_COUNTRY_AND_CANONICAL_NAME_CANDIDATE_NOT_GEO_VERIFIED"]
    z2020=[x for x in matches if x["Zenodo2020_actual_WPI_port_identity_verified"]]
    combined=[x for x in matches if x in cand and x in unct and x in z2020]
    audit={
        "schema":"chocho_2018_original_glsn_194_WPI_country_name_join_undirected_source_v1",
        "status":"FULL_2018_GLOBAL_LINER_EDGE_SOURCE_RECONSTRUCTED_AND_WPI194_CANDIDATES_ASSESSED",
        "original_GLSN_ports":len(ports),
        "original_GLSN_excel_edge_rows_directed":len(edges),
        "original_GLSN_unique_unordered_name_pairs":undirected,
        "original_GLSN_source_pairs_with_exactly_two_direction_rows":symmetric_name_pairs,
        "original_GLSN_pairs_with_wrong_mirror_or_weight":mismatched_reverse,
        "original_GLSN_pairs_with_missing_or_extra_direction_rows":asymmetric_counts,
        "original_GLSN_self_loop_rows":self_loops,
        "source_GLSN_edge_name_unknown_endpoints":dict(unknown_endpoints),
        "GLSN_source_countries_not_mapped_to_ISO2":dict(unresolved_country),
        "original_source_WPI_ports":194,"source_islands":51,
        "GLSN_unique_ISO3_ISO2_name_originalWPI_candidates":len(cand),
        "GLSN_candidate_original_island_groups":len({r["island_id"] for r in cand}),
        "GLSN_induced_modelled_undirected_WPI_candidate_edges":len(candidate_edges),
        "GLSN_original_directed_rows_in_candidate_edge_pairs":raw_directed_in_candidate,
        "GLSN_original_unordered_island_pairs_with_source_modelled_edges":len(modelled_between_island),
        "UNCTAD_name_country_port_candidates":len(unct),
        "UNCTAD_early_16of20_original_ports":sum(
            int(r["UNCTAD_2006_2010_available_quarters"])>=16 for r in unct),
        "original_AIS2020_verifiably_identified_WPI_ports":len(z2020),
        "candidate_ports_in_all_three_sources":len(combined),
        "all_three_source_original_WPI_IDs":[r["original_WPI_port_id"] for r in combined],
        "GLSN_has_NO_original_geographic_coords_and_WPI_cannot_be_verified_by_GLSN":True,
        "UNCTAD_provider_has_NO_geographic_coords_and_port_identity_not_geoverified":True,
        "GLSN_weight_NOT_actual_vessel_OD_or_colonization":True,
        "raw_source_XLSX_sha256":{"ports":source_meta["sources"][0]["original_sha256"],
                                 "edges":source_meta["sources"][1]["original_sha256"]},
        "out_sha256":{p.name:digest(p) for p in out.iterdir() if p.is_file()},
        "source_only_audit_time_utc":datetime.now(timezone.utc).isoformat(),
    }
    (out/"GLSN2018_WPI194_full_source_edge_audit.json").write_text(
        json.dumps(audit,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GLSN2018_WPI194_MODELLED_SOURCE_MATCH_JSON="+json.dumps({
        k:audit[k] for k in (
            "original_GLSN_ports","original_GLSN_excel_edge_rows_directed",
            "original_GLSN_unique_unordered_name_pairs",
            "original_GLSN_source_pairs_with_exactly_two_direction_rows",
            "original_GLSN_pairs_with_wrong_mirror_or_weight",
            "GLSN_unique_ISO3_ISO2_name_originalWPI_candidates",
            "GLSN_candidate_original_island_groups",
            "GLSN_induced_modelled_undirected_WPI_candidate_edges",
            "GLSN_original_unordered_island_pairs_with_source_modelled_edges",
            "UNCTAD_name_country_port_candidates",
            "candidate_ports_in_all_three_sources",
            "all_three_source_original_WPI_IDs"
        )},sort_keys=True),flush=True)
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    run(p.parse_args().out_dir)

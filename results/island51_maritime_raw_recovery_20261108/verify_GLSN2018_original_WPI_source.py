#!/usr/bin/env python3
"""Independent old-source XLSX/PLSCIQ/2020-AIS to WPI194 edge and island-pair audit.

Separate implementation from producing join: enumerates every original 2018
modelled source edge and reconstructs original WPI source-only pair support.
No biological observations or models are accessed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
import pycountry

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"results/island51_maritime_raw_recovery_20261108"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
ORIGINAL_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
UPSTREAM=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/original194_source_matched_port_coverage.csv"

def records(path):
    with Path(path).open("rt",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def normalize(label):
    raw=unicodedata.normalize("NFKD",str(label))
    tokens=re.findall("[A-Z0-9]+",raw.encode("ascii","ignore").decode("ascii").upper())
    if tokens[:2]==["PORT","OF"] and len(tokens)>2:
        tokens=tokens[2:]
    elif len(tokens)>1 and tokens[0]=="PORT":
        tokens=tokens[1:]
    while len(tokens)>1 and tokens[-1] in ("PORT","HARBOR","HARBOUR","TERMINAL"):
        tokens.pop()
    return "".join(tokens)

def iso(source):
    found=pycountry.countries.get(alpha_3=source.strip().upper())
    return found.alpha_2 if found else ""

def check():
    data=WPI.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if sha!=ORIGINAL_BLOB:
        raise ValueError("Source original WPI194 input changed")
    wpi=records(WPI)
    ports=records(SOURCE/"GLSN2018_ports_original_sheet_1.csv")
    edges=records(SOURCE/"GLSN2018_edges_original_sheet_1.csv")
    unctad=records(SOURCE/"UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv")
    ais=records(UPSTREAM)
    rec=json.loads((SOURCE/"GLSN2018_recovered_source_receipt.json").read_text())
    ure=json.loads((SOURCE/"UNCTAD_PLSCIQ_source_recovery_receipt.json").read_text())
    for part,filename in [("PORTS","GLSN2018_ports_original_sheet_1.csv"),("EDGES","GLSN2018_edges_original_sheet_1.csv")]:
        item=next(t for t in rec["sources"] if t["role"]==part)
        digest=hashlib.sha256((SOURCE/filename).read_bytes()).hexdigest()
        if len(item["sheets"])!=1 or item["sheets"][0]["extracted_csv_sha256"]!=digest:
            raise ValueError("Original GLSN2018 source SHA changed")
    if len(wpi)!=194 or len(ports)!=887 or len(edges)!=28104 or len(unctad)!=194 or len(ais)!=194:
        raise ValueError("Source port/network expected fixed 194/887/28104/194/194 rows")
    if not ure["source_version"]["is_exact_2026_10_30_original_recovered"]:
        raise ValueError("Original UNCTAD 2026-10-30 archived source vintage not exactly recovered")
    if ure["series"]!=1300 or ure["total_port_quarter_rows"]!=81900 or ure["nonmissing_port_quarter_rows"]!=57312:
        raise ValueError("Original UNCTAD historical index values no longer source-locked")
    count_glsn=Counter(normalize(p["port"]) for p in ports)
    provider=defaultdict(list)
    for row in ports:
        provider[(iso(row["ISO3"]),normalize(row["port"]))].append(row["port"])
    wpi_counts=Counter((p["country_code"],normalize(p["port_name"])) for p in wpi)
    expected={}
    name_to_wpi={}
    for w in wpi:
        pid=str(int(w["port_id"]))
        key=(w["country_code"],normalize(w["port_name"]))
        options=provider.get(key,[])
        yes=(len(options)==1 and count_glsn[normalize(options[0])]==1 and wpi_counts[key]==1)
        expected[pid]=yes
        if yes:
            if options[0] in name_to_wpi:raise ValueError("Two original WPI entries linked same provider name")
            name_to_wpi[options[0]]=pid
    byid={str(int(w["port_id"])):w for w in wpi}
    source_edge_undirected=defaultdict(list)
    for r in edges:
        a,b=r["port1"],r["port2"]
        if a==b:raise ValueError("Original GLSN self-loop unexpected")
        source_edge_undirected[tuple(sorted([a,b]))].append(r)
    # Always distinguish observed physical AIS voyages from 2018 liner co-service edges.
    matched={}
    island_edges=Counter()
    for key,rows in source_edge_undirected.items():
        a,b=key
        if a in name_to_wpi and b in name_to_wpi:
            i,j=name_to_wpi[a],name_to_wpi[b]
            if i==j:continue
            if (i,j) in matched or (j,i) in matched:
                raise ValueError("Duplicate source matched WPI unordered edge")
            matched[tuple(sorted((i,j),key=int))]=rows[0]
            ia,ib=byid[i]["island_id"],byid[j]["island_id"]
            if ia!=ib:
                island_edges[tuple(sorted((ia,ib),key=int))]+=1
    derived=records(Path(args.directory)/"GLSN2018_UNCTAD2006_2010_Zenodo2020_originalWPI194_candidate_join.csv")
    if len(derived)!=194 or {str(int(r["original_WPI_port_id"])) for r in derived}!=set(byid):
        raise ValueError("Incomplete original WPI194 derived cross-source panel")
    byu={str(int(r["original_WPI_id"])):r for r in unctad}
    bya={str(int(r["original_WPI_port_id"])):r for r in ais}
    for r in derived:
        pid=str(int(r["original_WPI_port_id"]))
        yes=expected[pid]
        actual=r["GLSN2018_status"]=="UNIQUE_ISO3_TO_ISO2_AND_NAME_CANDIDATE_NOT_GEO_VERIFIED"
        if yes!=actual or r["island_id"]!=byid[pid]["island_id"]:
            raise ValueError("GLSN 2018 original name/country source-only WPI matching changed")
        if r["GLSN2018_original_coords_INDEPENDENTLY_VERIFIED"].lower()!="false":
            raise ValueError("Original GLSN source has no lat/lon and may not be certified via geography")
        if (r["UNCTAD_PLSCIQ_match_status"]!=byu[pid]["status"] or
            (r["Zenodo2020_actual_WPI_port_identity_verified"].lower()=="true")!=(
                bya[pid]["verified_external_2020_endpoint"].lower()=="true")):
            raise ValueError("UNCTAD or AIS2020 original WPI joint comparison changed")
    out_edges=records(Path(args.directory)/"GLSN2018_originalWPI_source_candidate_model_edges.csv")
    if len(out_edges)!=len(matched):
        raise ValueError("Original GLSN modelled WPI candidate edge count differs from source")
    used=set()
    for r in out_edges:
        a,b=str(int(r["original_WPI_candidate_port_1"])),str(int(r["original_WPI_candidate_port_2"]))
        pair=tuple(sorted((a,b),key=int))
        if pair in used or pair not in matched:raise ValueError("Original WPI source modelled edge invalid/duplicated")
        used.add(pair)
        source=matched[pair]
        if (not math.isclose(float(r["original_source_weight_NOT_TRAFFIC_VOLUME"]),float(source["weight"]),rel_tol=1e-10) or
            not math.isclose(float(r["original_source_GC_km"]),float(source["GC"]),rel_tol=1e-10)):
            raise ValueError("Original GLSN edge weight/GC changed")
        if r["edge_is_liner_service_co_membership_not_actual_direct_voyage"].lower()!="true":
            raise ValueError("Co-service edge relabeled actual voyage")
    pairrows=records(Path(args.directory)/"GLSN2018_original51_all1275_modelled_island_pair_candidates.csv")
    if len(pairrows)!=1275:
        raise ValueError("All 1275 source island pair modelled candidate rows required")
    counted=set()
    for r in pairrows:
        ia,ib=r["original_island_a"],r["original_island_b"]
        pair=(ia,ib)
        if pair in counted or int(ia)>=int(ib):
            raise ValueError("Original WPI source island pair repeated/invalid")
        counted.add(pair)
        if int(r["GLSN2018_modelled_undirected_copresence_edges_in_matched_candidate_ports"])!=island_edges[pair]:
            raise ValueError("Original matched co-service source island edge count differs")
        if r["missing_model_source_edge_is_not_absence_of_ship_traffic"].lower()!="true":
            raise ValueError("Missing model link improperly asserted absence of vessel traffic")
    if len(counted)!=1275:
        raise ValueError("Original 51 landmass pairs not all checked")
    counts={
      "status":"EVERY_2018_EXTRACTED_SOURCE_EDGE_AND_ORIGINAL_WPI_1275_PAIR_RECONSTRUCTED",
      "original_GLSN_source_port_rows":len(ports),
      "original_GLSN_source_directed_rows":len(edges),
      "original_GLSN_source_unordered_edge_labels":len(source_edge_undirected),
      "source_GL_country_name_candidate_ports":sum(expected.values()),
      "source_GL_candidate_wpi_modelled_edges":len(matched),
      "source_GL_candidate_island_pairs":len(island_edges),
      "original_WPI_ports_with_other_sources_checked":len(derived),
      "original51_island_pairs_checked":len(counted),
      "original_UNCTAD_legacy_archive_exact_version":True,
      "no_true_shipping_or_butterfly_arrival_inference":True
    }
    (Path(args.directory)/"independent_full_GLSN2018_WPI_audit.json").write_text(json.dumps(counts,sort_keys=True,indent=2)+"\n")
    print("GLSN2018_INDEPENDENT_MODEL_SOURCE_AUDIT="+json.dumps(counts,sort_keys=True),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--directory",type=Path,required=True)
    args=ap.parse_args()
    check()

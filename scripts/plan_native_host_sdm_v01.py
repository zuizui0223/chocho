#!/usr/bin/env python3
"""Freeze native-screened host species and native-land climatic background, before LGM.

Preparation only: no 0BP model fit, LGM climatic score or genetic result opened.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

NATIVE_SOURCE_SHA="bf00cbbd49d901da80dcb30e3824c2ae30a427a17ace2b0e935b05f69fd72cd4"
HOST_MANIFEST_SHA="a28fba9efa8a9dc162336f62c883059c5dddbfabec2861f554d038dd98044".replace("162336","162336") # overwritten below by immutable expected digest
HOST_MANIFEST_SHA="a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
NATIVE_TAXA_SHA="8fbf5e0fab0df919862955f725a1822c7af05e3d59887d785cbc3d22b9b96dc0"
SUPPORT_SHA="d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f"
R=6371.0088
MAX_PRESENCES=200
MIN_PRESENCES=30
MIN_BACKGROUND=30
MAX_BACKGROUND=1000

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def rows(path):
    with Path(path).open("r",newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))

def stable_priority(p,tid):
    return hashlib.sha256(f"chocho-native-sdm-v01|{tid}|{p['gbif_id']}".encode()).hexdigest()

def spatial_thin(records,taxon,kilometers=10,cap=200):
    if not records:return []
    selected=[];xyz=[]
    theta=2*np.sin(kilometers/(2*R))
    for p in sorted(records,key=lambda row:(stable_priority(row,taxon),str(row["gbif_id"]))):
        lat=float(p["latitude"]);lon=float(p["longitude"])
        if not -90<=lat<=90 or not -180<=lon<=180:
            raise ValueError("bad WCVP-native source coordinate")
        phi=np.deg2rad(lat);lam=np.deg2rad(lon)
        v=np.array([np.cos(phi)*np.cos(lam),np.cos(phi)*np.sin(lam),np.sin(phi)])
        if xyz and np.min(np.linalg.norm(np.asarray(xyz)-v,axis=1))<theta-1e-12:
            continue
        selected.append(p);xyz.append(v)
        if len(selected)>=cap:break
    return selected

def make_plan(manifest,native_points,native_units,support):
    hosts=defaultdict(dict)
    for r in rows(manifest):
        key=str(int(float(r["wcvp_id_as_frozen"])))
        other=hosts[key].get("name")
        if other and other!=r["accepted_host_name"]:raise ValueError("conflicting WCVP accepted name")
        hosts[key]["name"]=r["accepted_host_name"]
    if len(hosts)!=518:raise ValueError("host taxon panel changed")
    native_codes=defaultdict(set)
    for r in rows(native_units):
        native_codes[str(int(r["accepted_plant_name_id"]))].add(r["area_code_l3"])
    if len(native_codes)!=517:raise ValueError("native WCVP taxon count changed")
    global_bg=rows(support)
    if len(global_bg)!=6373:raise ValueError("source WCVP representative land points changed")
    unit_points=defaultdict(list)
    for i,r in enumerate(global_bg):unit_points[r["area_code_l3"]].append(i)
    accepted=defaultdict(list)
    with gzip.open(native_points,"rt",newline="",encoding="utf8") as f:
        for p in csv.DictReader(f):
            if p["geography_status"]!="WITHIN_DOCUMENTED_NATIVE_L3":
                raise ValueError("unqualified WCVP-native observation in source")
            accepted[p["wcvp_id"]].append(p)
    if sum(map(len,accepted.values()))!=220137:raise ValueError("native observation source count changed")
    plants=[];selection=[]
    for taxon in sorted(hosts,key=lambda v:int(v)):
        chosen=spatial_thin(accepted.get(taxon,[]),taxon,kilometers=10,cap=MAX_PRESENCES)
        candidate_bg=sorted({i for code in native_codes.get(taxon,set()) for i in unit_points.get(code,[])})
        if len(candidate_bg)>MAX_BACKGROUND:
            candidate_bg=sorted(candidate_bg,key=lambda i:hashlib.sha256(
                f"chocho-native-background-v01|{taxon}|{i}".encode()).hexdigest())[:MAX_BACKGROUND]
            candidate_bg=sorted(candidate_bg)
        ready=len(chosen)>=MIN_PRESENCES and len(candidate_bg)>=MIN_BACKGROUND
        state=("READY_NATIVE_SDM" if ready else
               "NO_DOCUMENTED_NATIVE_WCVP" if taxon not in native_codes else
               "INSUFFICIENT_NATIVE_OBSERVATIONS" if len(chosen)<MIN_PRESENCES else
               "INSUFFICIENT_NATIVE_BACKGROUND")
        plants.append({"wcvp_id":taxon,"accepted_host_name":hosts[taxon]["name"],
            "native_region_L3_count":len(native_codes.get(taxon,set())),
            "original_native_observations":len(accepted.get(taxon,[])),
            "thinned_native_observations":len(chosen),
            "native_background_support_points":len(candidate_bg),
            "status":state})
        for p in chosen:
            selection.append({"wcvp_id":taxon,"gbif_id":p["gbif_id"],
                "latitude":p["latitude"],"longitude":p["longitude"],
                "native_wgsrpd3_code":p["native_wgsrpd3_code"],
                "uncertainty_m":p["uncertainty_m"]})
    census=dict(Counter(p["status"] for p in plants))
    return plants,selection,census

def main():
    ap=argparse.ArgumentParser()
    for p in ("host-manifest","native-points","native-units","support","out"):
        ap.add_argument("--"+p,required=True,type=Path)
    a=ap.parse_args()
    for path,expect in ((a.host_manifest,HOST_MANIFEST_SHA),
                        (a.native_points,NATIVE_SOURCE_SHA),
                        (a.native_units,NATIVE_TAXA_SHA),(a.support,SUPPORT_SHA)):
        if sha(path)!=expect:raise ValueError("source hash mismatch: "+str(path))
    plants,selected,counts=make_plan(a.host_manifest,a.native_points,a.native_units,a.support)
    a.out.mkdir(parents=True,exist_ok=True)
    with (a.out/"all_518_native_sdm_qualification.csv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(plants[0]),lineterminator="\n")
        w.writeheader();w.writerows(plants)
    with (a.out/"source_thinned_native_occurrences.csv.gz").open("wb") as raw:
        with gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0,compresslevel=6) as zipped:
            with io.TextIOWrapper(zipped,encoding="utf8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(selected[0]),lineterminator="\n")
                w.writeheader();w.writerows(selected)
    summary={"schema":"chocho_native_WCVP_host_SDM_0BP_LGM_eligibility_v0.1",
        "status":"SOURCE_FIXED_BEFORE_0BP_OR_LGM_MODELLING",
        "plant_host_taxa":518,"native_WCVP_taxa":517,
        "source_native_observation_points":220137,
        "thinned_plant_observations":len(selected),
        "eligibility_counts":counts,
        "model_eligible_taxa":counts.get("READY_NATIVE_SDM",0),
        "original_insect_species":641,"original_native_paired_insects":392,
        "source_sha256":{p.name:sha(p) for p in (a.host_manifest,a.native_points,a.native_units,a.support)},
        "model_rules":{"spatial_thinning_km":10,"max_plant_points":MAX_PRESENCES,
            "minimum_presences":MIN_PRESENCES,"minimum_native_background":MIN_BACKGROUND,
            "maximum_background_points":MAX_BACKGROUND,
            "background":"original 6373 WCVP support points inside each host's CURRENT documented native L3"},
        "no_paleoclimate_or_genetic_response_used":True,
        "limits":["Eligibility is based on current WCVP regional-native source and modern GBIF observations only",
          "Historical plant occupancy, dispersal and gene flow are not measured"]}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))

if __name__=="__main__":main()

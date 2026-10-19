#!/usr/bin/env python3
"""Qualify all original insects using native-WCVP-screened modern host references.

Present WCVP L3 ranges only: not historical plant distributions or genetic
results. Exact original 50km uncertainty-adjusted sampling exclusion, 10km
thinning, cap200; every original 641 species retained in the disposition table.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

NATIVE_FILE_SHA="bf00cbbd49d901da80dcb30e3824c2ae30a427a17ace2b0e935b05f69fd72cd4"
MANIFEST_SHA="a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
LOC_SHA="037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9"
EDGE_SHA="ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9"

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def load_csv(p):
    with Path(p).open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def qualifications(known, insect_n, points, native):
    source_ids=set(known)&set(native)
    observed_ids=set(p["wcvp_id"] for p in points)
    good=insect_n>=30 and len(points)>=30
    return {"paired":good,"all_native_hosts_observed":source_ids==set(known),
            "all_native_hosts_selected":observed_ids==set(known),
            "source_native_host_ids":len(source_ids),
            "selected_native_host_ids":len(observed_ids)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ttf-root",type=Path,required=True)
    p.add_argument("--native-observations",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    ttf=a.ttf_root
    sys.path[:0]=[str(ttf/"scripts"),str(ttf/"src")]
    from reference_helpers_min import join_hosts,filter_host_points,normalize_id
    from ttf.historical_host_memory import great_circle_quadrature_latlon
    hpath=ttf/"benchmarks/exploratory/lgm_external_reference_v01/host_reference_request_manifest.csv"
    ipath=ttf/"benchmarks/exploratory/lgm_external_reference_v01/all_641_species_external_reference.csv"
    original=ttf/"benchmarks/exploratory/lgm_yearstrat_398_climate_v03/all_641_eligibility_v0.3.csv"
    geopath=ttf/"input/native_geometry"
    lp=geopath/"historical_host_memory_localities_v0.1.csv"
    ep=geopath/"historical_host_memory_edges_v0.1.csv"
    for q,expected in [(hpath,MANIFEST_SHA),(lp,LOC_SHA),(ep,EDGE_SHA),
                       (a.native_observations,NATIVE_FILE_SHA)]:
        if sha(q)!=expected:raise ValueError("immutable source drift: "+str(q))
    refs={r["species"]:r for r in load_csv(ipath)}
    old={r["species"]:r for r in load_csv(original)}
    if not (len(refs)==len(old)==641):raise ValueError("original 641 cohort drift")
    hostsets=defaultdict(set)
    for r in load_csv(hpath):
        hostsets[r["insect_species"]].add(normalize_id(r["wcvp_id_as_frozen"]))
    if len(hostsets)!=641:raise ValueError("full 641 host identity drift")
    native=defaultdict(list)
    with gzip.open(a.native_observations,"rt",newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            if r["geography_status"]!="WITHIN_DOCUMENTED_NATIVE_L3":
                raise ValueError("native-only source has unqualified record")
            native[r["wcvp_id"]].append(r)
    if sum(map(len,native.values()))!=220137:raise ValueError("full native WCVP observation count drift")
    loc=defaultdict(dict);edge=defaultdict(list)
    for r in load_csv(lp):
        loc[r["species"]][int(r["locality_index"])]=(float(r["latitude"]),float(r["longitude"]))
    for r in load_csv(ep):
        edge[r["species"]].append((int(r["node_left"]),int(r["node_right"])))
    support={}
    for name,positions in loc.items():
        orig=np.array([positions[i] for i in sorted(positions)],float)
        pair=np.array([[positions[i],positions[j]] for i,j in edge[name]],float)
        interior=great_circle_quadrature_latlon(pair[:,0,:],pair[:,1,:],segment_points=5).reshape(-1,2)
        support[name]=np.vstack([orig,interior])
    if len(support)!=642:raise ValueError("original original-geography 642 entity mismatch")
    disposition=[];stored=[];c=Counter()
    for name in sorted(refs):
        known=hostsets[name]
        full=join_hosts(known & native.keys(),native)
        selected=filter_host_points(full,support[name],buffer_km=50,thinning_km=10,cap=200)
        insect_n=int(refs[name].get("external_reference_points_50km") or 0)
        q=qualifications(known,insect_n,selected,native)
        original_ok=old[name]["state"]=="READY_FOR_CLIMATE"
        transition=("retained" if q["paired"] and original_ok else
            "lost_native_screen" if original_ok and not q["paired"] else
            "new_native_only" if q["paired"] and not original_ok else "ineligible_both")
        c[transition]+=1
        c["insect_ge30"]+=int(insect_n>=30)
        c["native_host_ge30"]+=int(len(selected)>=30)
        c["paired_native_ge30"]+=int(q["paired"])
        c["paired_all_native_hosts_observed"]+=int(q["paired"] and q["all_native_hosts_observed"])
        c["paired_all_native_hosts_selected"]+=int(q["paired"] and q["all_native_hosts_selected"])
        disposition.append({"species":name,"insect_reference_n":insect_n,
            "native_host_quality_observations":len(full),
            "native_host_reference_n":len(selected),
            "known_host_taxa":len(known),"native_host_taxa_with_any_record":q["source_native_host_ids"],
            "native_host_taxa_in_thinned_cloud":q["selected_native_host_ids"],
            "paired_native_reference_ge30":q["paired"],
            "old_unfiltered_paired_ge30":original_ok,
            "all_known_hosts_observed_native":q["all_native_hosts_observed"],
            "all_known_hosts_retained_native":q["all_native_hosts_selected"],
            "change_from_original":transition})
        for item in selected:
            stored.append({"insect_species":name,"wcvp_id":item["wcvp_id"],
              "gbif_id":item["gbif_id"],"latitude":item["latitude"],
              "longitude":item["longitude"],"uncertainty_m":item["uncertainty_m"],
              "event_date":item.get("event_date",""),
              "native_wgsrpd3_code":item["native_wgsrpd3_code"],
              "nearest_evaluation_km":item["nearest_evaluation_km"]})
    if (c["retained"],c["lost_native_screen"],c["paired_native_ge30"])!=(392,6,392):
        raise ValueError("native source geometry cohort did not reproduce original")
    if (len(disposition),len(stored))!=(641,111759):
        raise ValueError("native selected reference counts differ from independently computed source")
    a.out.mkdir(parents=True,exist_ok=True)
    with (a.out/"all_641_native_reference_eligibility.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(disposition[0]),lineterminator="\n")
        w.writeheader();w.writerows(disposition)
    with (a.out/"all_641_native_reference_points.csv.gz").open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,filename="",mode="wb",mtime=0,compresslevel=6) as gz:
            with io.TextIOWrapper(gz,encoding="utf8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(stored[0]),lineterminator="\n")
                w.writeheader();w.writerows(stored)
    summary={"schema":"chocho_WCVP_native_screened_TTF_641_paired_reference_v01",
        "status":"MODERN_NATIVE_REGION_SCREEN_ONLY_NOT_LGM_OCCUPANCY",
        "original_species":641,"original_unfiltered_qualified":398,
        "native_qualified_species":c["paired_native_ge30"],
        "lost_after_native_screen":c["lost_native_screen"],
        "native_reference_point_rows":len(stored),
        "counts":dict(c),
        "source_sha256":{"native_points":NATIVE_FILE_SHA,
            "original_manifest":MANIFEST_SHA,"original_localities":LOC_SHA,
            "original_edges":EDGE_SHA},
        "genetic_response_used":False,
        "historical_plant_range_or_barrier_inferred":False,
        "limitations":["Current WCVP native L3 areas are broad, not 21ka occupancy",
           "Plant observations are nonrandom modern iNaturalist records",
           "Shared plant observations may occur across focal insect rows"]}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))

if __name__=="__main__":main()

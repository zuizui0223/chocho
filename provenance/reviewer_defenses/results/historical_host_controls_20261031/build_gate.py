#!/usr/bin/env python3
"""Recompute historical plant recording opportunity, NOT host absence or facilitation."""
from __future__ import annotations
from pathlib import Path
import csv, hashlib, json, io
ROOT=Path(__file__).resolve().parents[4]
B=Path("provenance/reviewer_defenses/results")
PATHS={"pairs":B/"island_host_network_20261030/island_host_network_pairs.csv",
       "plants":B/"island_dated_host_match_20261030/island_port_early_plant_availability.csv",
       "baseline":B/"maritime_host_steppingstone_20261030/near_island_record_order_baseline_20261030.json"}
PINS={"pairs":"5663d2ffe42419f5c4496b92c9f3bbf76175907b",
      "plants":"0a01831056c2cd59a4bd41bd1dbd4c84d77cf922",
      "baseline":"345c20c2c6492a29f5d23f06b5f95655f3e72f7c"}
ORDER=["PRIORITY_UNCLASSIFIED","FOCAL_HOST_OBSERVED","HIGH_RECORD_NO_NEW","LOW_BOTANICAL_RECORDS"]
FIELDS=["a","b","km","early_a","early_b","min_early","new_butterfly_pairs",
        "known_early_host_pairs","host_unknown_new_pairs","ordered_rerecord","gate"]
def inputs(root=ROOT):
    out={}
    for k,p in PATHS.items():
        data=(root/p).read_bytes()
        sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
        if sha!=PINS[k]:raise ValueError("Git source changed "+k)
        out[k]=data.decode("utf8")
    return out
def audit(root=ROOT):
    src=inputs(root)
    p=list(csv.DictReader(io.StringIO(src["pairs"])))
    plants=list(csv.DictReader(io.StringIO(src["plants"])))
    screen=json.loads(src["baseline"])
    if (len(p),len(plants),len(screen["pair_summaries"]))!=(1275,194,14):
        raise ValueError("Source population changed")
    sortedpair=lambda a,b:tuple(sorted([str(a),str(b)],key=int))
    early={}
    for x in plants:
        i=x["landmass_id"]
        early[i]=max(early.get(i,0),int(x["early_plant_records"]))
    base={sortedpair(x["island_a"],x["island_b"]):x for x in screen["pair_summaries"]}
    rows=[]
    for x in p:
        if float(x["port_distance_km"])>=200:continue
        a,b=sortedpair(x["island_a"],x["island_b"])
        s=base[a,b]
        if abs(float(x["port_distance_km"])-s["port_representative_distance_km"])>1e-9:
            raise ValueError("Island pair geography changed")
        n=int(x["shared_new_species"]);f=int(s["previously_prior_host_matched_focal_links"])
        if n!=int(s["shared_new_butterfly_species_pair_count"]) or f>n:
            raise ValueError("Butterfly host-pair count changed")
        ea,eb=early.get(a,0),early.get(b,0)
        m=min(ea,eb);other=n-f
        gate=("LOW_BOTANICAL_RECORDS" if m<100 else "FOCAL_HOST_OBSERVED" if f>0
              else "PRIORITY_UNCLASSIFIED" if other>0 else "HIGH_RECORD_NO_NEW")
        rows.append(dict(a=a,b=b,km=f'{float(x["port_distance_km"]):.6f}',
            early_a=ea,early_b=eb,min_early=m,new_butterfly_pairs=n,
            known_early_host_pairs=f,host_unknown_new_pairs=other,
            ordered_rerecord=int(s["upstream_re_record_before_downstream_first_year"]),
            gate=gate))
    rows.sort(key=lambda x:(ORDER.index(x["gate"]),-x["host_unknown_new_pairs"],
                            -x["min_early"],int(x["a"]),int(x["b"])))
    group={g:{"pairs":sum(x["gate"]==g for x in rows),
              "new":sum(x["new_butterfly_pairs"] for x in rows if x["gate"]==g),
              "unknown":sum(x["host_unknown_new_pairs"] for x in rows if x["gate"]==g)}
           for g in ORDER}
    if (len(rows),sum(x["new_butterfly_pairs"] for x in rows),
        group["PRIORITY_UNCLASSIFIED"]["pairs"],
        group["PRIORITY_UNCLASSIFIED"]["unknown"],
        group["FOCAL_HOST_OBSERVED"]["new"],group["LOW_BOTANICAL_RECORDS"]["new"]
        )!=(14,41,5,9,3,29):
        raise ValueError("Historical host comparison support changed")
    summary=dict(schema="chocho_historical_host_control_acquisition_gate_v1",
       status="EXPLORATORY_SELECTION_NOT_CAUSAL_COMPARISON",
       source_git_blob_pins=PINS,near_island_pairs=len(rows),
       unique_islands=len({i for x in rows for i in (x["a"],x["b"])}),
       shared_new_butterfly_pair_records=41,botanical_record_threshold_each_island_max_port=100,
       botanical_record_period="2001–2010",gate_counts=group,
       priority_pairs=[[x["a"],x["b"]] for x in rows if x["gate"]=="PRIORITY_UNCLASSIFIED"],
       no_confirmed_host_absence_controls=True,
       original_plant_species_dates_and_port_nonport_coordinates_missing=True,
       port_distances_are_2019_proxy_not_ship_routes=True)
    return rows,summary
def render_csv(rows):
    o=io.StringIO()
    w=csv.DictWriter(o,fieldnames=FIELDS,lineterminator="\n")
    w.writeheader();w.writerows(rows)
    return o.getvalue()
if __name__=="__main__":
    import argparse
    a=argparse.ArgumentParser()
    a.add_argument("--out-dir",type=Path,required=True)
    args=a.parse_args()
    rows,summary=audit()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"pair_priority.csv").write_text(render_csv(rows),encoding="utf8")
    (args.out_dir/"pair_priority_summary.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    print(summary["gate_counts"])

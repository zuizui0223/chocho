#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def load_panel(path):
    with path.open(newline="",encoding="utf-8") as f:
        names=sorted({str(r["species"]).strip() for r in csv.DictReader(f) if str(r.get("species") or "").strip()})
    if len(names)!=31: raise RuntimeError(f"expected 31 frozen occurrence species; got {len(names)}")
    return names

def load_pairs(path):
    by=defaultdict(set); meta={}
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sp=str(r["insect_species"]).strip(); pid=str(r["accepted_plant_name_id"]).strip()
            if not sp or not pid: continue
            by[sp].add(pid)
            meta[pid]={"accepted_name":str(r.get("accepted_name") or "").strip(),"family":str(r.get("family") or "").strip()}
    return by,meta

def load_units(path):
    out=defaultdict(set)
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            pid=str(r["accepted_plant_name_id"]).strip(); u=str(r["area_code_l3"]).strip()
            if pid and u: out[pid].add(u)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--occurrences-csv",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--native-distribution-csv",type=Path,required=True)
    ap.add_argument("--contemporary-distribution-csv",type=Path,required=True)
    ap.add_argument("--output-host-csv",type=Path,required=True)
    ap.add_argument("--output-cell-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()
    panel=load_panel(a.occurrences_csv); pairs,meta=load_pairs(a.insect_host_csv)
    native=load_units(a.native_distribution_csv); cont=load_units(a.contemporary_distribution_csv)
    host_rows=[]; cell_rows=[]; by_species={}
    for sp in panel:
        hosts=sorted(pairs.get(sp,set())); nu=set(); cu=set()
        for h in hosts: nu|=native.get(h,set()); cu|=cont.get(h,set())
        introduced_only=sorted(cu-nu); hrows=0
        for u in introduced_only:
            contributors=[h for h in hosts if u in cont.get(h,set()) and u not in native.get(h,set())]
            if not contributors: raise RuntimeError(f"no introduced contributor for {sp} {u}")
            hrows+=len(contributors)
            cell_rows.append({"species":sp,"wgsrpd3_code":u,"actual_host_count":len(contributors)})
            for h in contributors:
                host_rows.append({"species":sp,"wgsrpd3_code":u,"host_id":h,"host_name":meta[h]["accepted_name"],"host_family":meta[h]["family"]})
        by_species[sp]={"known_hosts":len(hosts),"native_resource_units":len(nu),"contemporary_resource_units":len(cu),"introduced_only_resource_units":len(introduced_only),"introduced_host_x_region_rows":hrows}
    if not host_rows or not cell_rows: raise RuntimeError("empty temporal resource panel")
    frozen_counts={
        "introduced_only_resource_cells":2035,
        "introduced_host_x_region_rows":5117,
        "unique_actual_hosts":215,
        "species_with_introduced_only_cells":29,
    }
    observed_counts={
        "introduced_only_resource_cells":len(cell_rows),
        "introduced_host_x_region_rows":len(host_rows),
        "unique_actual_hosts":len({r["host_id"] for r in host_rows}),
        "species_with_introduced_only_cells":sum(v["introduced_only_resource_units"]>0 for v in by_species.values()),
    }
    if observed_counts != frozen_counts:
        raise RuntimeError(f"frozen temporal resource panel drift: {observed_counts} != {frozen_counts}")
    a.output_host_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_host_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(host_rows[0])); w.writeheader(); w.writerows(host_rows)
    with a.output_cell_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(cell_rows[0])); w.writeheader(); w.writerows(cell_rows)
    counts=sorted(int(r["actual_host_count"]) for r in cell_rows)
    payload={"schema":"chocho_butterfly_resource_consumer_temporal_panel_v0.1","status":"OUTCOME_BLIND_INTRODUCED_ONLY_RESOURCE_PANEL","panel_species":len(panel),"introduced_only_resource_cells":len(cell_rows),"introduced_host_x_region_rows":len(host_rows),"unique_actual_hosts":len({r["host_id"] for r in host_rows}),"median_actual_hosts_per_cell":counts[len(counts)//2],"species_with_introduced_only_cells":sum(v["introduced_only_resource_units"]>0 for v in by_species.values()),"by_species":by_species,"definition":"All contemporary known-host resource regions outside each butterfly's native known-host union, using the frozen 31-species panel and frozen HOSTS-WCVP reconstruction. No butterfly future/outcome information is used.","claim_boundary":"This freezes the resource-side chronology panel only; no temporal result is computed here."}
    a.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in payload.items() if k!="by_species"},indent=2))
if __name__=="__main__": main()

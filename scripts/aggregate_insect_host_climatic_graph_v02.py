#!/usr/bin/env python3
"""Species-linked topology of potentially suitable 21ka native larval-host climate.

No observed ancient host occupancy, continuous fine-scale biological routes or
genetic responses are ever used. Does not rescue negative TTF host-identity null.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import defaultdict,Counter
from decimal import Decimal
from pathlib import Path
import numpy as np
from pmip4_potential_habitat_components_v02 import graph_metrics

ORIGINAL_HOST_MANIFEST="a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
SHAPE=(360,720)
BITS=360*720

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""):h.update(c)
    return h.hexdigest()

def accepted_id(value):
    n=Decimal(str(value))
    if not n.is_finite() or n<=0 or n!=n.to_integral_value():
        raise ValueError("bad source accepted plant ID")
    return str(int(n))

def host_model_status(host_ids,models):
    known=set(host_ids);modeled=set(models)
    if not known:return "NO_RECORDED_HOST"
    if known<=modeled:return "ALL_RECORDED_HOST_MODELS_AVAILABLE"
    if known&modeled:return "PARTIAL_RECORDED_HOST_MODELS_NOT_EVALUATED"
    return "NO_RECORDED_HOST_MODEL_AVAILABLE"

def union_known_host_potential(hosts,source_models,masked_now,masked_lgm):
    if not hosts or not set(hosts)<=set(source_models):
        raise ValueError("incomplete original host models cannot be silently treated as absent")
    order=[source_models[k] for k in sorted(set(hosts),key=lambda x:int(x))]
    now=np.any(np.unpackbits(np.asarray(masked_now)[order],axis=1)[:,:BITS].astype(bool),axis=0)
    then=np.any(np.unpackbits(np.asarray(masked_lgm)[order],axis=1)[:,:BITS].astype(bool),axis=0)
    return now.reshape(SHAPE),then.reshape(SHAPE)

def compare_component_metrics(mask0,mask21,lat):
    c0=graph_metrics(mask0,lat,4);c21=graph_metrics(mask21,lat,4)
    return {"area_km2_0BP":c0["suitable_km2"],"area_km2_21ka":c21["suitable_km2"],
        "area_km2_delta_21ka_minus_0BP":c21["suitable_km2"]-c0["suitable_km2"],
        "suitable_cells_0BP":c0["suitable_grid_cells"],
        "suitable_cells_21ka":c21["suitable_grid_cells"],
        "component_count_0BP":c0["components"],
        "component_count_21ka":c21["components"],
        "component_count_delta":c21["components"]-c0["components"],
        "largest_component_km2_0BP":c0["largest_component_km2"],
        "largest_component_km2_21ka":c21["largest_component_km2"],
        "largest_component_fraction_of_suitable_area_0BP":c0["fraction_area_in_largest_component"],
        "largest_component_fraction_of_suitable_area_21ka":c21["fraction_area_in_largest_component"],
        "sites_newly_suitable_21ka":int(np.count_nonzero(mask21&~mask0)),
        "sites_no_longer_suitable_21ka":int(np.count_nonzero(mask0&~mask21)),
        "sites_suitable_both":int(np.count_nonzero(mask0&mask21))}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--plant-mask",required=True,type=Path)
    p.add_argument("--plant-summary",required=True,type=Path)
    p.add_argument("--insect-host-manifest",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if sha(args.insect_host_manifest)!=ORIGINAL_HOST_MANIFEST:
        raise ValueError("source 641 insect WCVP associations changed")
    q=json.loads(args.plant_summary.read_text())
    if (q["status"]!="PHYSICAL_LAND_AND_CLIMATIC_POTENTIAL_GRAPH_NOT_LGM_OCCUPANCY_OR_GENETICS"
        or q["models_reproduced_on_original_6337_sites"]!=416
        or q["source_grid_cells"]!=BITS):
        raise ValueError("parent all-416 plant model graph not verified")
    links=defaultdict(set)
    with args.insect_host_manifest.open("r",newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            links[r["insect_species"]].add(accepted_id(r["wcvp_id_as_frozen"]))
    if len(links)!=641:raise ValueError("original 641 insect taxa lost")
    with np.load(args.plant_mask,allow_pickle=False) as z:
        ids=[accepted_id(k) for k in z["plant_wcvp_ids"]]
        m0=np.array(z["mask_0BP"],dtype="u1")
        m21=np.array(z["mask_21ka"],dtype="u1")
        shape=tuple(z["grid_shape"].tolist())
        lat=np.array(z["lat_centers"],float)
        if (len(ids)!=416 or len(set(ids))!=416 or shape!=SHAPE
            or m0.shape!=m21.shape!=(416,32400) or
            len(lat)!=360 or int(z["land_cut"])!=5 or
            abs(float(z["threshold_quantile"])-.1)>1e-9
            or not bool(z["nonanalog_excluded"])):
            raise ValueError("source models/masks/threshold not all source-bound")
        index={k:i for i,k in enumerate(ids)}
        model_ids=set(index)
        output=[];status=Counter()
        for species in sorted(links):
            known=links[species]
            s=host_model_status(known,model_ids)
            status[s]+=1
            row={"insect_species":species,"status":s,
                 "recorded_host_n":len(known),
                 "modelled_host_n":len(known&model_ids)}
            if s=="ALL_RECORDED_HOST_MODELS_AVAILABLE":
                now,old=union_known_host_potential(known,index,m0,m21)
                row.update(compare_component_metrics(now,old,lat))
            output.append(row)
    if status!=Counter({"ALL_RECORDED_HOST_MODELS_AVAILABLE":513,
                       "PARTIAL_RECORDED_HOST_MODELS_NOT_EVALUATED":83,
                       "NO_RECORDED_HOST_MODEL_AVAILABLE":45}):
        raise ValueError("original 513/83/45 host model qualification drift")
    complete=[x for x in output if x["status"]=="ALL_RECORDED_HOST_MODELS_AVAILABLE"]
    change=np.array([v["area_km2_delta_21ka_minus_0BP"] for v in complete],float)
    fragment=np.array([v["component_count_delta"] for v in complete],int)
    result={"schema":"chocho_513_insects_union_native_hosts_PMIP4_05deg_v02",
        "status":"HYPOTHETICAL_CLIMATIC_RESOURCE_CONNECTIVITY_NOT_REAL_LGM_OCCUPANCY",
        "original_insects":641,"all_recorded_hosts_modelled_insects":513,
        "source_native_plants":416,"statuses":dict(status),
        "area_decrease_species":int((change<0).sum()),
        "area_increase_species":int((change>0).sum()),
        "area_equal_species":int((change==0).sum()),
        "median_paired_21ka_minus_0BP_suitable_km2":float(np.median(change)),
        "component_count_increase_species":int((fragment>0).sum()),
        "component_count_decrease_species":int((fragment<0).sum()),
        "component_count_unchanged_species":int((fragment==0).sum()),
        "median_component_count_change":float(np.median(fragment)),
        "species_with_no_LGM_potential_climatic_habitat":int(sum(v["suitable_cells_21ka"]==0 for v in complete)),
        "source_sha256":{"plant_mask":sha(args.plant_mask),"plant_summary":sha(args.plant_summary),
                         "insect_original_hosts":ORIGINAL_HOST_MANIFEST},
        "limits":["The OR of all documented host's climate-suitable grid cells is not evidence of their actual occupancy or larval resource supply",
                  "4-neighbor components on coarse 0.5 degree PMIP4 physical land are potential substrate topology, not individual insect dispersal routes",
                  "Missing ancient climate analogues were censored; source sampling and current WCVP native ranges remain limiting",
                  "Model source's 1-3 documented larval hosts are not necessarily a complete physiological host catalogue",
                  "Multiple insects share host plant prediction surfaces and genetic response was never accessed"]}
    path=args.out/"all_641_insect_potential_21ka_host_components_v02.csv"
    cols=sorted(set().union(*(r.keys() for r in output)))
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,lineterminator="\n");w.writeheader();w.writerows(output)
    (args.out/"summary.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":main()

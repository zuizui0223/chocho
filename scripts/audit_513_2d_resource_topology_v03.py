#!/usr/bin/env python3
"""Source-locked robustness audit of 513 hypothetical 21ka climate-resource graphs.

Four/eight-neighbor sensitivity, area vs fragmentation, and dependence caused by
different insects using exactly the same recorded host-plant combination.
NO occupied LGM habitat, dispersal success or genetic response inferred.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path
import numpy as np
from pmip4_potential_habitat_components_v02 import graph_metrics

EXPECTED={
    "plant_masks":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
    "original_insects":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044",
    "original513":"27cd6d45ec67382d7606b40f74e232ac3bccbd3e7d0118498635b58d3de27cad"
}
GRID=(360,720)
N=360*720

def sha256_file(path):
    sha=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda:fh.read(1024*1024),b""):
            sha.update(chunk)
    return sha.hexdigest()

def accepted(value):
    d=Decimal(str(value))
    if not d.is_finite() or d<=0 or d!=d.to_integral_value():
        raise ValueError("invalid accepted host ID")
    return str(int(d))

def load_csv(path):
    with Path(path).open(newline="",encoding="utf8") as file:
        return list(csv.DictReader(file))

def summarize(rows):
    if not rows:raise ValueError("no climate-resource graph rows")
    def freq(test):return int(sum(bool(test(row)) for row in rows))
    def med(key):return float(np.median([float(row[key]) for row in rows]))
    return {"n":len(rows),
        "potential_suitable_area_declines":freq(lambda x:x["area_delta_km2"]<0),
        "four_neighbor_components_decline":freq(lambda x:x["component_delta_4"]<0),
        "eight_neighbor_components_decline":freq(lambda x:x["component_delta_8"]<0),
        "four_neighbor_largest_fraction_declines":freq(lambda x:x["largest_fraction_delta_4"]<0),
        "eight_neighbor_largest_fraction_declines":freq(lambda x:x["largest_fraction_delta_8"]<0),
        "four_neighbor_largest_area_declines":freq(lambda x:x["largest_area_delta_4_km2"]<0),
        "eight_neighbor_largest_area_declines":freq(lambda x:x["largest_area_delta_8_km2"]<0),
        "four_neighbor_components_and_largest_fraction_both_decline":freq(
           lambda x:x["component_delta_4"]<0 and x["largest_fraction_delta_4"]<0),
        "eight_neighbor_components_and_largest_fraction_both_decline":freq(
           lambda x:x["component_delta_8"]<0 and x["largest_fraction_delta_8"]<0),
        "paired_median_component_delta_4":med("component_delta_4"),
        "paired_median_component_delta_8":med("component_delta_8"),
        "paired_median_largest_fraction_delta_4":med("largest_fraction_delta_4"),
        "paired_median_largest_fraction_delta_8":med("largest_fraction_delta_8"),
        "paired_median_largest_area_delta_4_km2":med("largest_area_delta_4_km2")}

def compute_masks(plant_indices,masked0,masked21):
    ind=np.asarray(plant_indices,int)
    if len(ind)==0 or len(set(ind))!=len(ind):raise ValueError("incomplete/duplicated host indices")
    def unpack(matrix):
        mask=np.unpackbits(matrix[ind],axis=1)
        return np.any(mask[:,:N].astype(bool),axis=0).reshape(GRID)
    return unpack(masked0),unpack(masked21)

def group_metrics(mask0,mask21,lat):
    z={}
    for neighborhood in (4,8):
        a=graph_metrics(mask0,lat,connectivity=neighborhood)
        b=graph_metrics(mask21,lat,connectivity=neighborhood)
        assert a["fraction_area_in_largest_component"] is not None
        assert b["fraction_area_in_largest_component"] is not None
        z.update({
            f"component_count_0BP_{neighborhood}":a["components"],
            f"component_count_21ka_{neighborhood}":b["components"],
            f"component_delta_{neighborhood}":b["components"]-a["components"],
            f"largest_fraction_0BP_{neighborhood}":a["fraction_area_in_largest_component"],
            f"largest_fraction_21ka_{neighborhood}":b["fraction_area_in_largest_component"],
            f"largest_fraction_delta_{neighborhood}":b["fraction_area_in_largest_component"]-a["fraction_area_in_largest_component"],
            f"largest_area_delta_{neighborhood}_km2":b["largest_component_km2"]-a["largest_component_km2"]})
        if neighborhood==4:
            z["area_delta_km2"]=b["suitable_km2"]-a["suitable_km2"]
            z["area_0BP_km2"]=a["suitable_km2"]
            z["area_21ka_km2"]=b["suitable_km2"]
    return z

def audit(plant_file,original_insects,original513):
    for key,path in (("plant_masks",plant_file),("original_insects",original_insects),
                     ("original513",original513)):
        if sha256_file(path)!=EXPECTED[key]:
            raise ValueError("source SHA mismatch "+key)
    original=load_csv(original513)
    if len(original)!=641 or len(set(r["insect_species"] for r in original))!=641:
        raise ValueError("original all-641 dispositions changed")
    by_sp={r["insect_species"]:r for r in original}
    links=defaultdict(set)
    for x in load_csv(original_insects):
        links[x["insect_species"]].add(accepted(x["wcvp_id_as_frozen"]))
    if set(links)!=set(by_sp):raise ValueError("source insect cohort mismatch")
    qualified={s for s,r in by_sp.items() if r["status"]=="ALL_RECORDED_HOST_MODELS_AVAILABLE"}
    if len(qualified)!=513:raise ValueError("not 513 originally qualified insect species")
    combos=defaultdict(list)
    for s in qualified:
        combos[tuple(sorted(links[s],key=int))].append(s)
    if len(combos)!=349:raise ValueError("non-independent hostset count not recovered")
    rows=[]
    with np.load(plant_file,allow_pickle=False) as data:
        plants=[accepted(x) for x in data["plant_wcvp_ids"]]
        if len(plants)!=len(set(plants))!=416:raise ValueError("original plant model cohort drift")
        ref={key:i for i,key in enumerate(plants)}
        m0,m21=data["mask_0BP"],data["mask_21ka"]
        lat=np.asarray(data["lat_centers"],float)
        if not (m0.shape==m21.shape==(416,32400) and
                tuple(data["grid_shape"].tolist())==GRID and
                abs(float(data["threshold_quantile"])-.10)<1e-10 and
                int(data["land_cut"])==5 and bool(data["nonanalog_excluded"])):
            raise ValueError("source masks are not original prespecified primary data")
        for count,(hosts,species) in enumerate(sorted(combos.items())):
            if not all(host in ref for host in hosts):
                raise ValueError("originally qualifying insect lacks a known source host")
            masks=compute_masks([ref[x] for x in hosts],m0,m21)
            measured=group_metrics(*masks,lat)
            for sp in sorted(species):
                old=by_sp[sp]
                for e in ("0BP","21ka"):
                    if int(old["component_count_"+e])!=measured["component_count_"+e+"_4"]:
                        raise ValueError("original 4-neighbor connected component count did not replay: "+sp)
                    if abs(float(old["area_km2_"+e])-measured["area_"+e+"_km2"])>1e-4:
                        raise ValueError("original 4-neighbor suitable area did not replay: "+sp)
                    if abs(float(old["largest_component_fraction_of_suitable_area_"+e])
                           -measured["largest_fraction_"+e+"_4"])>1e-7:
                        raise ValueError("original largest cluster fraction did not replay: "+sp)
                rows.append({"insect_species":sp,"recorded_host_n":len(hosts),
                    "hostset_key":";".join(hosts),"hostset_insect_species_n":len(species),**measured})
            if (count+1)%75==0:
                print(json.dumps({"distinct_hostsets_verified":count+1,"total":349}),flush=True)
    if len(rows)!=513:raise ValueError("lost insect results")
    rows.sort(key=lambda x:x["insect_species"])
    groups={r["hostset_key"]:r for r in rows}
    unique=list(groups.values())
    report={"schema":"chocho_native_host_513_insects_349_hostsets_4vs8_neighbor_audit_v03",
        "status":"EXPLORATORY_POTENTIAL_CLIMATE_TOPOLOGY_SENSITIVITY_NOT_GENETICS",
        "original_insects":641,"full_host_model_insects":513,
        "distinct_recorded_host_combinations":len(unique),
        "distinct_host_plant_models_used":len(set().union(*(set(combo) for combo in combos))),
        "shared_hosts_cause_nonindependent_insect_results":True,
        "recorded_host_count_distribution":dict(Counter(r["recorded_host_n"] for r in rows)),
        "source4_original_insect_values_verified":True,
        "species_weighted":summarize(rows),
        "unique_hostset_equal_weight":summarize(unique),
        "suitable_area_contracting_species":summarize([r for r in rows if r["area_delta_km2"]<0]),
        "suitable_area_expanding_species":summarize([r for r in rows if r["area_delta_km2"]>0]),
        "source_sha256":dict(EXPECTED),
        "limitations":[
            "All 513 insects mapped to 349 different recorded host combinations; no independent 513 species resource graph replicates",
            "Eight neighbors permit diagonal corners and do not prove migratory plant corridors",
            "Shrinking suitable components can reduce number of components without increasing fraction concentrated in largest remaining core",
            "Modern native-region trained climate niche fits and exclusion of past nonanalog climates are substantial untested assumptions",
            "Potentially passable 0.5deg terrain is not actual ancient host range or taxon-specific dispersal",
            "No independent LGM insect genetic response opened; previous true-host-vs-other-plant negative screening preserved"]}
    return rows,report

def main():
    p=argparse.ArgumentParser()
    for field in ("plant-mask","insect-host-manifest","original513","out"):
        p.add_argument("--"+field,type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    rows,report=audit(args.plant_mask,args.insect_host_manifest,args.original513)
    with (args.out/"all_513_4v8_topology_and_hostset.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    (args.out/"summary.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()

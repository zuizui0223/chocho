#!/usr/bin/env python3
"""Apply verified PMIP4 physical land/ice masks to 416 plant and 641 insect scores.

The result is a climate-suitable AND ice-free *source-sampled site* comparison,
NOT plant occupancy, continuous climate corridors or a genetic barrier effect.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import defaultdict,Counter
from decimal import Decimal
from pathlib import Path
import numpy as np

SOURCES={
"plant_model":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733",
"plant_npz":"8559ec905cddd1a0b9579530551c9980ab3519daa92902df150a04687bf81159",
"insect_host":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def rows(path):
    with Path(path).open("r",encoding="utf8",newline="") as stream:
        return list(csv.DictReader(stream))

def accepted_id(value):
    d=Decimal(str(value).strip())
    if not d.is_finite() or d<=0 or d!=d.to_integral_value():raise ValueError("invalid plant ID")
    return str(int(d))

def any_recorded_host_opportunity(now,old,thresholds,free_now,free_old):
    now=np.asarray(now,float);old=np.asarray(old,float)
    thresholds=np.asarray(thresholds,float)
    free_now=np.asarray(free_now,bool);free_old=np.asarray(free_old,bool)
    if (now.ndim!=2 or old.shape!=now.shape or thresholds.shape!=(len(now),) or
        free_now.shape!=free_old.shape!=(now.shape[1],) or
        not np.isfinite(now).all() or not np.isfinite(old).all() or
        not np.isfinite(thresholds).all()):
        raise ValueError("source model, 10th percentile threshold, or physical site mask invalid")
    clim0=np.any(now>=thresholds[:,None],axis=0)
    clim21=np.any(old>=thresholds[:,None],axis=0)
    phys0=clim0&free_now;phys21=clim21&free_old
    return {
       "climatic_suitable_fraction_0ka":float(np.mean(clim0)),
       "climatic_suitable_fraction_21ka":float(np.mean(clim21)),
       "climatic_delta_21ka_minus_0ka":float(np.mean(clim21)-np.mean(clim0)),
       "icefree_climatic_fraction_0ka":float(np.mean(phys0)),
       "icefree_climatic_fraction_21ka":float(np.mean(phys21)),
       "icefree_climatic_delta_21ka_minus_0ka":float(np.mean(phys21)-np.mean(phys0)),
       "climate_suitable_but_physically_blocked_0ka_n":int(np.count_nonzero(clim0&~free_now)),
       "climate_suitable_but_physically_blocked_21ka_n":int(np.count_nonzero(clim21&~free_old)),
       "icefree_suitable_site_gain_21ka_n":int(np.count_nonzero(~phys0&phys21)),
       "icefree_suitable_site_loss_21ka_n":int(np.count_nonzero(phys0&~phys21))}

def delta_stats(rows,field):
    d=np.array([float(r[field]) for r in rows])
    if not np.isfinite(d).all():raise ValueError("nonfinite source-specific site proportion")
    return {"species":int(len(d)),"decrease":int((d<0).sum()),"increase":int((d>0).sum()),
            "unchanged":int((d==0).sum()),"median_paired_delta":float(np.median(d))}

def main():
    parser=argparse.ArgumentParser()
    for name in ("plant-model","plant-npz","host-associations","land-ice-mask","source-physics-summary","out"):
        parser.add_argument("--"+name,required=True,type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for key,path in (("plant_model",args.plant_model),("plant_npz",args.plant_npz),
                     ("insect_host",args.host_associations)):
        if sha(path)!=SOURCES[key]:raise ValueError("native-WCVP plant or original insect-host source digest drift")
    source=json.loads(args.source_physics_summary.read_text())
    if (source["status"]!="PHYSICAL_LAND_ICE_TOPOLOGY_NOT_21KA_PLANT_HABITAT_OR_GENETICS" or
        source["source_sha256"]["21ka"]!="3decf87e11a2d2afe65a7358d904214aea2e5c59ec92ab8227f36070beb98e94" or
        source["support_ice_free_land_0ka"]!=6047 or source["support_ice_free_land_21ka"]!=5391):
        raise ValueError("original source-verified PMIP4 physical masks changed")
    plants=rows(args.plant_model)
    if len(plants)!=518:raise ValueError("original 518 accepted plants changed")
    with np.load(args.plant_npz,allow_pickle=False) as score, np.load(args.land_ice_mask,allow_pickle=False) as pmip:
        ids=[accepted_id(x) for x in score["plant_wcvp_id"]]
        if ids!=[accepted_id(p["wcvp_id"]) for p in plants]:
            raise ValueError("plant accepted identity/score index changed")
        now=score["scores_present"];glacial=score["scores_21ka"]
        if now.shape!=glacial.shape!=(518,6337) or not np.array_equal(
            score["global_support_indices"],pmip["global_support_indices"]):
            raise ValueError("source point or model-score coordinate correspondence drift")
        free0=pmip["source_age0_ice_free"];free21=pmip["source_age21_ice_free"]
        if free0.shape!=free21.shape!=(6337,) or not np.isin(free0,[0,1]).all() or not np.isin(free21,[0,1]).all():
            raise ValueError("PMIP4 physical passability source does not contain binary 6337 points")
        model={id_:i for i,(id_,p) in enumerate(zip(ids,plants))
               if p["model_status"]=="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"}
        if len(model)!=416:raise ValueError("native plant model cohort qualification changed")
        threshold=np.array([float(p["threshold_0BP_presence_quantile_10"])
                            if id_ in model else np.nan for id_,p in zip(ids,plants)])
        plant_results=[]
        for id_,i in sorted(model.items(),key=lambda x:int(x[0])):
            plant_results.append({"wcvp_id":id_,"accepted_name":plants[i]["accepted_host_name"],
               **any_recorded_host_opportunity(now[[i]],glacial[[i]],threshold[[i]],free0,free21)})
        host_links=defaultdict(set)
        for r in rows(args.host_associations):
            sp=r["insect_species"].strip()
            if not sp:raise ValueError("unidentified insect")
            host_links[sp].add(accepted_id(r["wcvp_id_as_frozen"]))
        if len(host_links)!=641:raise ValueError("all original 641 insects not represented")
        results=[];statuses=Counter()
        for sp,hostset in sorted(host_links.items()):
            if hostset <= model.keys():status="COMPLETE"
            elif hostset & model.keys():status="PARTIAL_UNSCORED"
            else:status="NONE_UNSCORED"
            statuses[status]+=1
            row={"insect_species":sp,"status":status,
                 "recorded_host_n":len(hostset),"available_model_n":len(hostset&model.keys())}
            if status=="COMPLETE":
                inds=[model[h] for h in sorted(hostset,key=int)]
                row.update(any_recorded_host_opportunity(now[inds],glacial[inds],threshold[inds],free0,free21))
            results.append(row)
    if dict(statuses)!={"COMPLETE":513,"PARTIAL_UNSCORED":83,"NONE_UNSCORED":45}:
        raise ValueError("fixed 513 insects' complete host model availability changed")
    plant_old=delta_stats(plant_results,"climatic_delta_21ka_minus_0ka")
    plant_new=delta_stats(plant_results,"icefree_climatic_delta_21ka_minus_0ka")
    computed=[r for r in results if r["status"]=="COMPLETE"]
    insect_old=delta_stats(computed,"climatic_delta_21ka_minus_0ka")
    insect_new=delta_stats(computed,"icefree_climatic_delta_21ka_minus_0ka")
    for name,values in (("all_416_plants_physically_masked_v01.csv",plant_results),
                        ("all_641_insects_physically_masked_v01.csv",results)):
        columns=sorted(set().union(*(r.keys() for r in values)))
        with (args.out/name).open("w",newline="",encoding="utf8") as stream:
            writer=csv.DictWriter(stream,fieldnames=columns,lineterminator="\n")
            writer.writeheader();writer.writerows(values)
    q={"schema":"chocho_PMIP4_glacial_terrain_gated_climatic_opportunity_v01",
       "status":"POTENTIAL_CLIMATIC_SITES_WITH_PHYSICAL_LAND_ICE_GATING_NOT_CONNECTIVITY",
       "current_native_support_sites":6337,"source_plant_models":416,
       "original_insect_panel":641,"insects_with_all_recorded_host_models":513,
       "insect_availability_status_counts":dict(statuses),
       "plants_without_geography_mask":plant_old,"plants_with_geography_mask":plant_new,
       "insects_without_geography_mask":insect_old,"insects_with_geography_mask":insect_new,
       "no_true_historical_host_occupancy_or_genetic_response_used":True,
       "source_hashes":dict(SOURCES),
       "limitations":["Site masks use 6337 WCVP points designed from current native geography and do not cover newly exposed 21ka shelves",
         "True climatic-corridor passability needs a continuous 2D climate-habitat and dispersal model, not isolated support sites",
         "The original 0BP-fitted niche models often extrapolate to nonanalog 21ka conditions",
         "No host-species-specific genetic divergence test, ancient occupancy or migration inferred",
         "This does not reopen the frozen negative alternative-plant TTF v0.5 result"]}
    (args.out/"summary.json").write_text(json.dumps(q,indent=2,sort_keys=True)+"\n")
    print(json.dumps(q,sort_keys=True))
if __name__=="__main__":main()

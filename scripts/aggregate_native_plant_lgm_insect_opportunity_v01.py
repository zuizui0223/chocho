#!/usr/bin/env python3
"""Map fitted 0BP/21ka host climate suitability to potential insect resource sites.

This is descriptive macroecological *climatic opportunity*; it is not occupied
21ka host habitat, dispersal connectivity, or observed insect genetic response.
The binary OR across recorded hosts is not a calibrated probability of host use.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path
import numpy as np

EXPECTED={
    "model_summary":"c345e52f893aa23a2ed6ccb0a675b56db92f1c9f246cb662c866cbce85b916c0",
    "plant_scores":"8559ec905cddd1a0b9579530551c9980ab3519daa92902df150a04687bf81159",
    "plant_models_csv":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733",
    "original_insect_host":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"
}

def sha256(p:Path)->str:
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def idstr(value):
    d=Decimal(str(value).strip())
    if not d.is_finite() or d<=0 or d!=d.to_integral_value():
        raise ValueError("bad frozen WCVP accepted plant identity")
    return str(int(d))

def load_rows(path):
    with Path(path).open("r",newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def classify_insect_host_models(host_ids,ready):
    # Any unavailable accepted host invalidates a complete-resource inference.
    k=set(host_ids)
    if not k:return "NO_RECORDED_HOSTS"
    if k<=ready:return "COMPLETE_KNOWN_HOST_MODELS"
    if k&ready:return "PARTIAL_KNOWN_HOST_MODELS_NOT_SCORED"
    return "NO_KNOWN_HOST_MODELS"

def evaluate_any_host_climate(site_now,site_lgm,thresholds):
    p=np.asarray(site_now,float);q=np.asarray(site_lgm,float)
    t=np.asarray(thresholds,float)
    if (p.ndim!=2 or q.shape!=p.shape or len(p)!=len(t)
        or len(p)==0 or not np.isfinite(p).all()
        or not np.isfinite(q).all() or not np.isfinite(t).all()):
        raise ValueError("nonfinite or nonparallel host climate reference surface")
    now=np.any(p>=t[:,None],axis=0)
    lgm=np.any(q>=t[:,None],axis=0)
    n=now.size;ever=np.count_nonzero(now|lgm);both=np.count_nonzero(now&lgm)
    return {"sites":n,
            "present_suitable_sites":int(now.sum()),
            "LGM_suitable_sites":int(lgm.sum()),
            "retained_both_sites":int(both),
            "lost_since_LGM_sites":int(np.count_nonzero(now&~lgm)),
            "gain_in_LGM_sites":int(np.count_nonzero(lgm&~now)),
            "present_suitable_fraction":float(now.mean()),
            "LGM_suitable_fraction":float(lgm.mean()),
            "delta_LGM_minus_present":float(lgm.mean()-now.mean()),
            "climatic_site_Jaccard":float(both/ever) if ever else None}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for field in ("plant-summary","plant-models","plant-scores","insect-hosts","out"):
        ap.add_argument("--"+field,required=True,type=Path)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    pinned={"model_summary":args.plant_summary,"plant_models_csv":args.plant_models,
            "plant_scores":args.plant_scores,"original_insect_host":args.insect_hosts}
    for name,path in pinned.items():
        if sha256(path)!=EXPECTED[name]:
            raise ValueError(f"original source hash mismatch: {name}")
    summary=json.loads(args.plant_summary.read_text())
    if (summary["status"]!="EXPLORATORY_CLIMATIC_SUITABILITY_NOT_ACTUAL_LGM_PLANT_OCCUPANCY"
        or summary["model_complete_taxa"]!=416 or
        summary["global_sites_comparable_0BP_vs_21ka"]!=6337
        or summary["no_lgm_occupancy_or_genetic_response_claim"] is not True):
        raise ValueError("native 416-plant model result source drift")
    plant_rows=load_rows(args.plant_models)
    if len(plant_rows)!=518 or len({r["wcvp_id"] for r in plant_rows})!=518:
        raise ValueError("all 518 plant records not in original order")
    with np.load(args.plant_scores,allow_pickle=False) as q:
        ids=[idstr(x) for x in q["plant_wcvp_id"]]
        now=q["scores_present"]
        then=q["scores_21ka"]
        sites=q["global_support_indices"]
        if not (len(ids)==len(plant_rows)==518 and
                ids==[idstr(p["wcvp_id"]) for p in plant_rows] and
                now.shape==then.shape==(518,6337) and len(sites)==6337 and
                np.all(np.diff(sites)>0)):
            raise ValueError("native plant row/score/site matrix alignment mismatch")
        id_to_ix={value:i for i,value in enumerate(ids)}
        if len(id_to_ix)!=518:raise ValueError("duplicate model identity")
        modelled={r["wcvp_id"] for r in plant_rows
                  if r["model_status"]=="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"}
        if len(modelled)!=416:raise ValueError("416 source-qualified model outcomes changed")
        thresholds=np.array([float(r["threshold_0BP_presence_quantile_10"])
                             if r["model_status"]=="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"
                             else np.nan for r in plant_rows],float)
        if not all(np.isfinite(thresholds[id_to_ix[k]]) for k in modelled):
            raise ValueError("missing current fitted native threshold")
        if not np.all(np.isfinite(now[[id_to_ix[k] for k in modelled]])) or not np.all(np.isfinite(then[[id_to_ix[k] for k in modelled]])):
            raise ValueError("model-classified plant climate scores unexpectedly nonfinite")
        links=defaultdict(set)
        for r in load_rows(args.insect_hosts):
            sp=r["insect_species"].strip()
            host=idstr(r["wcvp_id_as_frozen"])
            if not sp or host not in id_to_ix:
                raise ValueError("unknown original insect or plant host")
            links[sp].add(host)
        if len(links)!=641:
            raise ValueError("641 insects with original host links not recovered")
        result=[];c=Counter()
        for species in sorted(links):
            hostset=links[species]
            state=classify_insect_host_models(hostset,modelled)
            c[state]+=1
            record={"insect_species":species,"status":state,
                    "recorded_host_taxa_n":len(hostset),
                    "modelled_host_taxa_n":len(hostset&modelled),
                    "unmodelled_host_taxa_n":len(hostset-modelled)}
            if state=="COMPLETE_KNOWN_HOST_MODELS":
                ids_for_sp=sorted(hostset,key=lambda k:int(k))
                idx=[id_to_ix[k] for k in ids_for_sp]
                record.update(evaluate_any_host_climate(now[idx],then[idx],thresholds[idx]))
            result.append(record)
    if (c["COMPLETE_KNOWN_HOST_MODELS"],
        c["PARTIAL_KNOWN_HOST_MODELS_NOT_SCORED"],
        c["NO_KNOWN_HOST_MODELS"])!=(513,83,45):
        raise ValueError("non-response 513/83/45 insect-host model availability changed")
    valid=[r for r in result if r["status"]=="COMPLETE_KNOWN_HOST_MODELS"]
    neg=sum(r["delta_LGM_minus_present"]<0 for r in valid)
    pos=sum(r["delta_LGM_minus_present"]>0 for r in valid)
    # No fixed biological-signal sign acceptance threshold: recompute descriptively.
    delta=np.array([r["delta_LGM_minus_present"] for r in valid])
    summary_out={"schema":"chocho_native_host_climate_potential_per_insect_v0.1",
         "status":"EXPLORATORY_CLIMATIC_RESOURCE_OPPORTUNITY_NOT_LGM_OCCUPANCY_OR_GENETICS",
         "original_insects":641,"recorded_plant_hosts":518,
         "native_host_climate_models_complete":416,"common_climate_support_sites":6337,
         "source_suitability_threshold":"each plant model 0BP native-presence 10th percentile",
         "species_status_counts":dict(c),"fully_modelled_insects":len(valid),
         "species_LGM_climate_opportunity_decrease":neg,
         "species_LGM_climate_opportunity_increase":pos,
         "species_LGM_climate_opportunity_equal":len(valid)-neg-pos,
         "median_species_delta_LGM_minus_present":float(np.median(delta)),
         "median_species_present_suitable_fraction":float(np.median([r["present_suitable_fraction"] for r in valid])),
         "median_species_LGM_suitable_fraction":float(np.median([r["LGM_suitable_fraction"] for r in valid])),
         "species_share_negative":float(neg/len(valid)),
         "source_sha256":dict(EXPECTED),
         "limitations":["Union of host climatic suitability masks is not plant presence, insect presence, larval resources actually available or resource graph connectivity",
           "Source host lists include only 1 to 3 accepted taxa and are not independently validated complete larval host ranges",
           "All 6337 reference sites are equal weighted point locations, not area-weighted global pixels",
           "Many historical plant climates are nonanalog and modern WCVP native L3 is coarse",
           "Insects share host plants and belong to dependent genera, so species counts are not 513 independent model fits",
           "No genetic response, glacial ice mask, sea-level connectivity or actual historical plant occurrence inferred"]}
    path=args.out/"all_641_insect_climatic_resource_opportunity.csv"
    keys=sorted(set().union(*(r.keys() for r in result)))
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader();w.writerows(result)
    (args.out/"summary.json").write_text(json.dumps(summary_out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary_out,sort_keys=True))

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Reconstruct the FROZEN 0BP-fitted native-plant logistic function from source 0BP scores.

Model restoration uses ONLY the original 0BP predicted probabilities and its
corresponding original 0BP CHELSA environmental design matrix. Archived 21ka
probabilities are held out for an exact 1e-5 model-replay validation. Novel
21ka climate locations are NOT known historical occupancy or genetic evidence.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.special import expit
from sklearn.preprocessing import PolynomialFeatures
from run_native_plant_21ka_climate_hindcast_v01 import read_csv,climate_sample
from fit_native_plant_05deg_glacial_potential_graph_v02 import read_grid
from pmip4_potential_habitat_components_v02 import graph_metrics

VARS=("bio01","bio07","bio12","bio15")
SHA={
    "models":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733",
    "scores":"8559ec905cddd1a0b9579530551c9980ab3519daa92902df150a04687bf81159",
    "main_mask":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
    "support":"d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f",
    "pmip05":"d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db",
}
GRID=(360,720)
N=360*720
REPLAY_TOL=1e-5

def digest(path):
    sha=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):sha.update(chunk)
    return sha.hexdigest()

def polynomial_features(x):
    arr=np.asarray(x,float)
    if arr.ndim!=2 or arr.shape[1]!=4 or not np.isfinite(arr).all():
        raise ValueError("invalid source 4-var CHELSA predictor")
    poly=PolynomialFeatures(2,include_bias=True)
    features=poly.fit_transform(arr)
    if features.shape[1]!=15:raise ValueError("quadratic source model feature count drift")
    return features

def frozen_logit_coefficients(source_features,present_scores):
    """Reverse original 15-parameter logistic function; 21ka is not used to fit."""
    features=np.asarray(source_features,float)
    probs=np.asarray(present_scores,float)
    if (features.ndim!=2 or features.shape[1]!=15 or len(features)!=len(probs)
        or not np.isfinite(features).all() or not np.isfinite(probs).all()
        or np.any((probs<0)|(probs>1))):
        raise ValueError("source 0BP model probabilities invalid")
    keep=(probs>1e-5)&(probs<1-1e-5)
    if keep.sum()<60 or np.linalg.matrix_rank(features[keep])!=15:
        raise ValueError("source frozen present model not identifiable in available 0BP scores")
    p=probs[keep]
    y=np.log(p)-np.log1p(-p)
    weight=np.sqrt(p*(1-p))
    c=np.linalg.lstsq(features[keep]*weight[:,None],y*weight,rcond=1e-12)[0]
    if c.shape!=(15,) or not np.isfinite(c).all():
        raise ValueError("non-finite source model parameters")
    return c,int(keep.sum())

def probabilities(features,coefficients):
    return expit(np.asarray(features,float)@np.asarray(coefficients,float))

def summarize_pair(main0,main21,extended0,extended21,lat):
    if (main0.shape!=main21.shape!=GRID or extended0.shape!=extended21.shape!=GRID):
        raise ValueError("unmatched source climate masks")
    original=[graph_metrics(main0,lat,4),graph_metrics(main21,lat,4)]
    extended=[graph_metrics(extended0,lat,4),graph_metrics(extended21,lat,4)]
    def report(x):
        a,b=x
        return {"area_delta_km2":b["suitable_km2"]-a["suitable_km2"],
            "component_delta":b["components"]-a["components"],
            "largest_fraction_delta":(
                b["fraction_area_in_largest_component"]-
                a["fraction_area_in_largest_component"]
                if a["fraction_area_in_largest_component"] is not None and
                   b["fraction_area_in_largest_component"] is not None else None),
            "largest_area_delta_km2":b["largest_component_km2"]-a["largest_component_km2"]}
    return {"main":report(original),"extrapolated":report(extended)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for field in ("source-models","source-scores","source-primary-mask","global-support",
                  "pmip05","chelsa","out"):
        parser.add_argument("--"+field,required=True,type=Path)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for key,path in [("models",a.source_models),("scores",a.source_scores),
                     ("main_mask",a.source_primary_mask),("support",a.global_support),
                     ("pmip05",a.pmip05)]:
        if digest(path)!=SHA[key]:raise ValueError("original model/source identity changed: "+key)
    models=read_csv(a.source_models)
    support=read_csv(a.global_support)
    if len(models)!=518 or len(support)!=6373:
        raise ValueError("original WCVP host model or reference geography incomplete")
    with np.load(a.source_scores,allow_pickle=False) as archive:
        original_ids=[str(k) for k in archive["plant_wcvp_id"]]
        support_ids=np.array(archive["global_support_indices"],int)
        saved0=np.asarray(archive["scores_present"],float)
        saved21=np.asarray(archive["scores_21ka"],float)
        if (original_ids!=[r["wcvp_id"] for r in models] or
            saved0.shape!=saved21.shape!=(518,6337) or
            len(support_ids)!=6337):
            raise ValueError("original 518 WCVP model predictions do not align")
    with np.load(a.source_primary_mask,allow_pickle=False) as archive:
        original_mask_ids=[str(k) for k in archive["plant_wcvp_ids"]]
        original_mask0=np.asarray(archive["mask_0BP"],dtype="u1")
        original_mask21=np.asarray(archive["mask_21ka"],dtype="u1")
        lat=np.asarray(archive["lat_centers"],float)
        if (original_mask0.shape!=original_mask21.shape!=(416,32400)
            or len(original_mask_ids)!=416 or not bool(archive["nonanalog_excluded"])
            or tuple(archive["grid_shape"].tolist())!=GRID):
            raise ValueError("original main climate-analog-only masked landscape changed")
    with np.load(a.pmip05,allow_pickle=False) as pmip:
        free0=np.asarray(pmip["pmip0_icefree_subcells"],dtype="u1").reshape(-1)
        free21=np.asarray(pmip["pmip21_icefree_subcells"],dtype="u1").reshape(-1)
        finite_old0=np.asarray(pmip["CHELSA0_fourvar_finite"],bool).reshape(-1)
        finite_old21=np.asarray(pmip["CHELSA21_fourvar_finite"],bool).reshape(-1)
        if not (free0.shape==free21.shape==finite_old0.shape==finite_old21.shape==(N,)):
            raise ValueError("PMIP4 original 0.5deg climate geometry drift")
    c0=[a.chelsa/f"CHELSA_TraCE21K_{v}_20_V1.0.tif" for v in VARS]
    c21=[a.chelsa/f"CHELSA_TraCE21K_{v}_-190_V1.0.tif" for v in VARS]
    if any(not x.is_file() for x in c0+c21):raise FileNotFoundError("8 original CHELSA source rasters needed")
    globalcoords=[(r["longitude"],r["latitude"]) for r in support]
    orig_current=climate_sample(c0,globalcoords)
    orig_lgm=climate_sample(c21,globalcoords)
    good=np.isfinite(orig_current).all(axis=1)
    comparable=good & np.isfinite(orig_lgm).all(axis=1)
    if good.sum()!=6343 or not np.array_equal(np.flatnonzero(comparable),support_ids):
        raise ValueError("source 6343/6337 reference-site climate identity changed")
    mu=orig_current[good].mean(axis=0);sigma=orig_current[good].std(axis=0)
    if not np.all(sigma>0):raise ValueError("source climate standardization undefined")
    source_F0=polynomial_features((orig_current[support_ids]-mu)/sigma)
    source_F21=polynomial_features((orig_lgm[support_ids]-mu)/sigma)
    grid0=read_grid(c0).reshape(N,4);grid21=read_grid(c21).reshape(N,4)
    finite0=np.isfinite(grid0).all(axis=1);finite21=np.isfinite(grid21).all(axis=1)
    if not np.array_equal(finite0,finite_old0) or not np.array_equal(finite21,finite_old21):
        raise ValueError("source 0.5 degree CHELSA coverage mismatch")
    active0=finite0&(free0>=5);active21=finite21&(free21>=5)
    F0=polynomial_features((grid0[active0]-mu)/sigma)
    F21=polynomial_features((grid21[active21]-mu)/sigma)
    rows=[];new0=[];new21=[];idx={v:i for i,v in enumerate(original_mask_ids)}
    successful=0;max_error0=0.;max_error21=0.
    for j,plant in enumerate(models):
        tid=original_ids[j]
        entry={"wcvp_id":tid,"accepted_name":plant["accepted_host_name"],
               "original_status":plant["model_status"]}
        if plant["model_status"]!="COMPLETE_POTENTIAL_CLIMATE_HINDCAST":
            entry["status"]="ORIGINAL_PLANT_MODEL_NOT_AVAILABLE"
            rows.append(entry);continue
        if tid not in idx:raise ValueError("plant originally modeled but no v0.2 baseline mask")
        coeff,n_train=frozen_logit_coefficients(source_F0,saved0[j])
        error0=float(np.max(np.abs(probabilities(source_F0,coeff)-saved0[j])))
        error21=float(np.max(np.abs(probabilities(source_F21,coeff)-saved21[j])))
        max_error0=max(max_error0,error0);max_error21=max(max_error21,error21)
        if error0>REPLAY_TOL or error21>REPLAY_TOL:
            diagnostic={"plant":tid,"name":plant["accepted_host_name"],
                "source_0bp_error":error0,"source_21ka_holdout_error":error21,
                "source_present_n":n_train}
            print(json.dumps({"SOURCE_LOGIT_REPLAY_FAILED":diagnostic},sort_keys=True),flush=True)
            raise ValueError("original 0BP model and 21ka heldout scores cannot be reproduced")
        threshold=float(plant["threshold_0BP_presence_quantile_10"])
        if not 0<threshold<1:raise ValueError("bad original frozen suitability threshold")
        model0=np.zeros(N,bool);model21=np.zeros(N,bool)
        model0[active0]=probabilities(F0,coeff)>=threshold
        model21[active21]=probabilities(F21,coeff)>=threshold
        k=idx[tid]
        original0=np.unpackbits(original_mask0[k])[:N].astype(bool)
        original21=np.unpackbits(original_mask21[k])[:N].astype(bool)
        missing_from_replay0=int(np.count_nonzero(original0&~model0))
        missing_from_replay21=int(np.count_nonzero(original21&~model21))
        # Superset of frozen v0.2 climate-analog-only model mask: never erase
        # a previously supported original location through inverse-logit rounding.
        extended0=model0|original0
        extended21=model21|original21
        assert not np.any(original0&~extended0) and not np.any(original21&~extended21)
        old2d0=original0.reshape(GRID);old2d21=original21.reshape(GRID)
        ext2d0=extended0.reshape(GRID);ext2d21=extended21.reshape(GRID)
        metrics=summarize_pair(old2d0,old2d21,ext2d0,ext2d21,lat)
        entry.update(status="VERIFIED_FROZEN_LOGIT_EXTRAPOLATION_SENSITIVITY",
                     source_0bp_score_max_error=error0,
                     heldout_21ka_score_max_error=error21,
                     source_0bp_midrange_sites=n_train,
                     source_primary_cells_not_reproduced_near_threshold_0bp=missing_from_replay0,
                     source_primary_cells_not_reproduced_near_threshold_21ka=missing_from_replay21,
                     cells_added_to_primary_mask_0bp=int(np.count_nonzero(extended0&~original0)),
                     cells_added_to_primary_mask_21ka=int(np.count_nonzero(extended21&~original21)),
                     original_climate_analog_only=metrics["main"],
                     extrapolation_included=metrics["extrapolated"])
        rows.append(entry)
        new0.append((tid,np.packbits(extended0)))
        new21.append((tid,np.packbits(extended21)))
        successful+=1
        if successful%60==0:
            print(json.dumps({"source_0BP_logistic_models_restored":successful,
                              "max_0bp_error":max_error0,
                              "max_21ka_holdout_error":max_error21}),flush=True)
    if successful!=416 or len(rows)!=518 or [v[0] for v in new0]!=original_mask_ids:
        raise ValueError("not all original 416 models or WCVP identities reconstituted")
    np.savez_compressed(a.out/"plant_potential_climate_masks_extrapolation_included_v03.npz",
          plant_wcvp_ids=np.asarray(original_mask_ids,dtype="<U16"),
          mask_0BP=np.stack([r[1] for r in new0]),
          mask_21ka=np.stack([r[1] for r in new21]),
          grid_shape=np.asarray(GRID),lat_centers=lat.astype("float32"),
          land_cut=np.array(5),threshold_quantile=np.array(.10),
          nonanalog_excluded=np.array(False),
          original_source_main_mask_sha256=np.array(SHA["main_mask"]))
    with (a.out/"all_518_native_host_model_replay_and_extrapolation.csv").open("w",newline="",encoding="utf8") as file:
        cols=sorted(set().union(*(r.keys() for r in rows)))
        writer=csv.DictWriter(file,fieldnames=cols,lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,dict) else v
                             for k,v in r.items()})
    modeled=[r for r in rows if r["status"]=="VERIFIED_FROZEN_LOGIT_EXTRAPOLATION_SENSITIVITY"]
    main=np.array([r["original_climate_analog_only"]["area_delta_km2"] for r in modeled])
    extra=np.array([r["extrapolation_included"]["area_delta_km2"] for r in modeled])
    mfrac=[r["original_climate_analog_only"]["largest_fraction_delta"] for r in modeled]
    efrac=[r["extrapolation_included"]["largest_fraction_delta"] for r in modeled]
    report={"schema":"chocho_source_0BP_logit_reconstruction_416_native_host_05deg_extrapolation_v03",
       "status":"SOURCE_VERIFIED_EXTRAPOLATION_UNCERTAINTY_NOT_OCCUPIED_LGM_HABITAT",
       "modeled_plant_taxa":416,"input_plant_taxa":518,
       "source_replay_max_absolute_error_0bp":max_error0,
       "source_holdout_max_absolute_error_21ka":max_error21,
       "source_21ka_predictions_used_only_for_validation":True,
       "reconstructed_source_model":"quadratic4var logistic logit, 15 terms; coefficients inferred from 0BP saved source probabilities; 21ka independently verified",
       "original_0BP_logistic_fit_data_used":False,
       "original_416_main_masks_inherited_byte_for_byte":True,
       "original_main_plant_area_decrease":int(np.sum(main<0)),
       "extrapolation_included_plant_area_decrease":int(np.sum(extra<0)),
       "original_main_largest_component_fraction_decrease":int(sum(x is not None and x<0 for x in mfrac)),
       "extrapolation_included_largest_component_fraction_decrease":int(sum(x is not None and x<0 for x in efrac)),
       "original_main_median_paired_area_delta_km2":float(np.median(main)),
       "extrapolation_included_median_paired_area_delta_km2":float(np.median(extra)),
       "source_main_threshold_border_disagree_0BP":int(sum(r["source_primary_cells_not_reproduced_near_threshold_0bp"] for r in modeled)),
       "source_main_threshold_border_disagree_21ka":int(sum(r["source_primary_cells_not_reproduced_near_threshold_21ka"] for r in modeled)),
       "source_sha256":dict(SHA),
       "interpretation":["Extrapolation scenario only: predictions in finite source 21ka climates outside observed modern environmental support are not validated occupancy",
         "The original 416 main masks are preserved and never refitted or retuned",
         "Inverse-logit coefficients fit ONLY saved 0BP probabilities. Original saved 21ka probabilities are external model identity checks",
         "0.5-degree physical land/ice approximation, nonrandom plant sampling and unknown ancient dispersal constrain all conclusions",
         "No genetic response, historical observed plant occurrence, biological interaction effect or migration inferred"]}
    (a.out/"summary.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()

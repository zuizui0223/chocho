#!/usr/bin/env python3
"""Reproduce 416 frozen native-WCVP 0BP SDMs and project potential habitat to 0.5°.

All outputs are source-pinned *climatic* suitability under PMIP4 ice-free land
constraints, not real 21 ka plant occupancy, dispersal, historical refugia or
insect genetic barriers. Inputs and eligibility are fixed before projection.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures

from pmip4_potential_habitat_components_v02 import graph_metrics
from run_native_plant_21ka_climate_hindcast_v01 import (
    read_csv,climate_sample,background_support_indices
)

ROOT_SHA={
 "plan":"1b009db9cc813a5ab992dd1baace23ea1aeed3673dbae5064a642be876c489c0",
 "native_points":"f580dfc76c7eb3d2b7ffa68e54e35bf6534578465d40e997f5bcfec0f19efc79",
 "native_l3":"8fbf5e0fab0df919862955f725a1822c7af05e3d59887d785cbc3d22b9b96dc0",
 "support":"d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f",
 "original_models":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733",
 "original_scores":"8559ec905cddd1a0b9579530551c9980ab3519daa92902df150a04687bf81159",
 "coarse_pmip_climate":"d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db",
 "coarse_pmip_summary":"2c18201c522906e7872617488d7382ad5d64ddb57741237024bd7fb488c1a15d",
}
VARS=("bio01","bio07","bio12","bio15")
THRESHOLDS=(.05,.10,.15)
GRID_SHAPE=(360,720)
GRID_N=360*720
SCORE_TOL=1e-5

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as inp:
        for chunk in iter(lambda:inp.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def read_grid(assets):
    transform=from_origin(-180,90,.5,.5)
    z=[]
    for p in assets:
        with rasterio.open(p) as src:
            if not src.crs or not src.crs.is_geographic:
                raise ValueError("original CHELSA source not geographic")
            with WarpedVRT(src,crs="EPSG:4326",transform=transform,
                           width=720,height=360,resampling=Resampling.nearest) as vrt:
                field=vrt.read(1,masked=True).astype(np.float32)
                z.append(np.asarray(field.filled(np.nan),dtype=np.float32))
    result=np.stack(z,axis=-1)
    if result.shape!=(360,720,4):raise ValueError("climate 0.5 degree grid misaligned")
    return result

def fit_original_classifier(present,background):
    p=np.asarray(present,float);b=np.asarray(background,float)
    if p.ndim!=2 or b.ndim!=2 or p.shape[1:]!=(4,) or b.shape[1:]!=(4,) or len(p)<30 or len(b)<30:
        raise ValueError("original native plant/background training insufficient")
    if not np.isfinite(p).all() or not np.isfinite(b).all():
        raise ValueError("nonfinite native plant source")
    train=np.vstack([p,b]);response=np.r_[np.ones(len(p)),np.zeros(len(b))]
    polynomial=PolynomialFeatures(degree=2,include_bias=False)
    x=polynomial.fit_transform(train)
    if x.shape[1]!=14:raise ValueError("original four climate variables quadratic degree mismatch")
    clf=LogisticRegression(C=.3,solver="liblinear",class_weight="balanced",
                           max_iter=1000,random_state=0,tol=1e-5)
    clf.fit(x,response)
    return clf,polynomial

def model_probs(model,poly,climate):
    a=np.asarray(climate,float)
    if a.ndim!=2 or a.shape[1]!=4 or not np.isfinite(a).all():
        raise ValueError("nonfinite projected climate")
    return model.predict_proba(poly.transform(a))[:,1]

def grid_validity_and_analog(climate,allowed,training):
    """Do not turn missing/nonanalog LGM cells into absence or suitable habitat."""
    a=np.asarray(climate,float);mask=np.asarray(allowed,bool)
    training=np.asarray(training,float)
    if a.shape!=(GRID_N,4) or mask.shape!=(GRID_N,) or training.ndim!=2 or training.shape[1]!=4:
        raise ValueError("incompatible grid climate vs physical source")
    available=np.isfinite(a).all(axis=1)&mask
    lo=training.min(axis=0);hi=training.max(axis=0)
    pad=.05*np.maximum(hi-lo,.001)
    analog=np.zeros(GRID_N,bool)
    x=a[available]
    analog[available]=np.all((x>=lo-pad)&(x<=hi+pad),axis=1)
    return available,analog

def suitable_mask(probs,threshold,free_subcells,analog,climate_finite,land_cut=5,include_nonanalog=False):
    probs=np.asarray(probs,float)
    if (probs.shape!=(GRID_N,) or np.shape(free_subcells)!=(GRID_N,) or
        np.shape(analog)!=(GRID_N,) or np.shape(climate_finite)!=(GRID_N,)):
        raise ValueError("0.5 degree terrain, climate, and plant projection are not aligned")
    allowed=(np.asarray(free_subcells)>=land_cut)&climate_finite
    if not include_nonanalog:allowed &= analog
    return allowed & np.isfinite(probs) & (probs>=float(threshold))

def count_available(physical,finite,analog):
    mask=np.asarray(physical,bool);finite=np.asarray(finite,bool);analog=np.asarray(analog,bool)
    return {"icefree_cells":int(mask.sum()),"icefree_with_climate_n":int((mask&finite).sum()),
        "icefree_without_climate_n":int((mask&~finite).sum()),
        "icefree_climate_nonanalog_n":int((mask&finite&~analog).sum()),
        "icefree_climate_analog_n":int((mask&finite&analog).sum())}

def main():
    p=argparse.ArgumentParser()
    for key in ("plant-plan","thinned-native","native-units","support","source-models",
                "source-scores","pmip-coarse","pmip-summary","chelsa-root","out"):
        p.add_argument("--"+key,required=True,type=Path)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for key,path in (("plan",args.plant_plan),("native_points",args.thinned_native),
                     ("native_l3",args.native_units),("support",args.support),
                     ("original_models",args.source_models),("original_scores",args.source_scores),
                     ("coarse_pmip_climate",args.pmip_coarse),("coarse_pmip_summary",args.pmip_summary)):
        if sha(path)!=ROOT_SHA[key]:raise ValueError("source changed after 0.5 degree graph contract: "+key)
    source_summary=json.loads(args.pmip_summary.read_text())
    if source_summary["status"]!="SOURCE_FIRST_PALEOLAND_CLIMATE_COVERAGE_NOT_SPECIES_CONNECTIVITY":
        raise ValueError("PMIP paleoland climate preflight not source complete")
    plants=read_csv(args.plant_plan);old_models=read_csv(args.source_models)
    if len(plants)!=len(old_models)!=518:raise ValueError("518 source host model identities changed")
    if [r["wcvp_id"] for r in plants]!=[r["wcvp_id"] for r in old_models]:
        raise ValueError("original native plant order changed")
    support=read_csv(args.support)
    if len(support)!=6373:raise ValueError("original source background coordinates missing")
    with np.load(args.pmip_coarse,allow_pickle=False) as v:
        free0=np.asarray(v["pmip0_icefree_subcells"],np.uint8).reshape(-1)
        free21=np.asarray(v["pmip21_icefree_subcells"],np.uint8).reshape(-1)
        preflight0=np.asarray(v["CHELSA0_fourvar_finite"],bool).reshape(-1)
        preflight21=np.asarray(v["CHELSA21_fourvar_finite"],bool).reshape(-1)
    if not all(a.shape==(GRID_N,) for a in (free0,free21,preflight0,preflight21)):
        raise ValueError("source compatibility grid shape")
    present_paths=[args.chelsa_root/f"CHELSA_TraCE21K_{v}_20_V1.0.tif" for v in VARS]
    ancient_paths=[args.chelsa_root/f"CHELSA_TraCE21K_{v}_-190_V1.0.tif" for v in VARS]
    if any(not x.is_file() for x in present_paths+ancient_paths):
        raise FileNotFoundError("source 0BP/21ka CHELSA missing")
    grid0=read_grid(present_paths).reshape(GRID_N,4)
    grid21=read_grid(ancient_paths).reshape(GRID_N,4)
    finite0=np.isfinite(grid0).all(axis=1);finite21=np.isfinite(grid21).all(axis=1)
    if not np.array_equal(preflight0,finite0) or not np.array_equal(preflight21,finite21):
        raise ValueError("model input CHELSA 0.5 degree grid disagrees with SOURCE FIRST preflight")
    pres=defaultdict(list)
    with gzip.open(args.thinned_native,"rt",newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):pres[r["wcvp_id"]].append(r)
    if sum(len(v) for v in pres.values())!=73386:raise ValueError("source native plant presence N drift")
    region=defaultdict(set)
    for r in read_csv(args.native_units):
        region[str(int(r["accepted_plant_name_id"]))].add(r["area_code_l3"])
    if len(region)!=517:raise ValueError("WCVP species native region drift")
    coords=[];lookup={}
    def ix(items):
        indices=[]
        for lon,lat in items:
            key=(round(float(lon),8),round(float(lat),8))
            if key not in lookup:lookup[key]=len(coords);coords.append(key)
            indices.append(lookup[key])
        return np.asarray(indices,int)
    global_idx=ix([(r["longitude"],r["latitude"]) for r in support])
    plant_index={}
    for plant in plants:
        key=plant["wcvp_id"]
        if plant["status"]=="READY_NATIVE_SDM":
            plant_index[key]=ix([(r["longitude"],r["latitude"]) for r in pres[key]])
    if len(plant_index)!=416:raise ValueError("qualified 416 plant source N drift")
    climate=climate_sample(present_paths,coords)
    previous=climate_sample(ancient_paths,[coords[i] for i in global_idx])
    global_source=climate[global_idx]
    valid0=np.isfinite(global_source).all(axis=1)
    if int(valid0.sum())!=6343:raise ValueError("0BP support-climate mean SD source mismatch")
    mu=global_source[valid0].mean(axis=0);sd=global_source[valid0].std(axis=0)
    if not np.all(sd>0):raise ValueError("climate standard deviation invalid")
    global_now=(global_source-mu)/sd
    global_then=(previous-mu)/sd
    comparable=valid0 & np.isfinite(previous).all(axis=1)
    if int(comparable.sum())!=6337:raise ValueError("original 6337 comparable source reference site N drift")
    reference_now=global_now[comparable];reference_then=global_then[comparable]
    with np.load(args.source_scores,allow_pickle=False) as old:
        original_plant_ids=[str(k) for k in old["plant_wcvp_id"]]
        old_idx=old["global_support_indices"].astype(int)
        source_scores_now=old["scores_present"]
        source_scores_then=old["scores_21ka"]
        if (original_plant_ids!=[p["wcvp_id"] for p in plants] or
            not np.array_equal(old_idx,np.flatnonzero(comparable)) or
            source_scores_now.shape!=source_scores_then.shape!=(518,6337)):
            raise ValueError("original 416 plant model reconstruction source alignment mismatch")
        # Retain arrays in process memory until every species is verified.
        source_scores_now=np.array(source_scores_now)
        source_scores_then=np.array(source_scores_then)
    scaled0=(grid0-mu)/sd;scaled21=(grid21-mu)/sd
    # Use original physical source 1/9 subcell passability as superset for all checks.
    relevant0=finite0&(free0>=1);relevant21=finite21&(free21>=1)
    relevant0_idx=np.flatnonzero(relevant0);relevant21_idx=np.flatnonzero(relevant21)
    lat=89.75-.5*np.arange(360)
    main_masks0=[];main_masks21=[]
    result=[];max_source_error=0.;max_threshold_error=0.
    status=Counter()
    for j,plant in enumerate(plants):
        key=plant["wcvp_id"]
        record={"wcvp_id":key,"accepted_host_name":plant["accepted_host_name"],
                "source_native_region_n":int(plant["native_region_L3_count"]),
                "source_native_presence_n":int(plant["thinned_native_observations"]),
                "source_native_background_n":int(plant["native_background_support_points"])}
        if plant["status"]!="READY_NATIVE_SDM":
            record["status"]=plant["status"];result.append(record);status[record["status"]]+=1;continue
        pb=plant_index[key]
        pi=np.asarray((climate[pb]-mu)/sd,float)
        bg_ids=background_support_indices(key,region,support)
        if len(bg_ids)!=int(plant["native_background_support_points"]):
            raise ValueError("frozen WCVP background source selection drift")
        bi=global_now[bg_ids]
        fp=np.isfinite(pi).all(axis=1);fb=np.isfinite(bi).all(axis=1)
        if fp.sum()<30 or fb.sum()<30:raise ValueError("formerly successful native training now inadequate")
        data_p=pi[fp];data_b=bi[fb]
        clf,poly=fit_original_classifier(data_p,data_b)
        predicted_original_now=model_probs(clf,poly,reference_now)
        predicted_original_then=model_probs(clf,poly,reference_then)
        if not np.isfinite(source_scores_now[j]).all() or not np.isfinite(source_scores_then[j]).all():
            raise ValueError("original successfully modeled species has nonfinite verification score")
        err=max(float(np.max(np.abs(predicted_original_now-source_scores_now[j]))),
                float(np.max(np.abs(predicted_original_then-source_scores_then[j]))))
        max_source_error=max(max_source_error,err)
        if err>SCORE_TOL:raise ValueError(f"original 6337 source SDM prediction not replayed: {key} {err}")
        fitted_presence=model_probs(clf,poly,data_p)
        original_threshold=float(old_models[j]["threshold_0BP_presence_quantile_10"])
        threshold=float(np.quantile(fitted_presence,.10))
        max_threshold_error=max(max_threshold_error,abs(threshold-original_threshold))
        if abs(threshold-original_threshold)>SCORE_TOL:
            raise ValueError(f"original source threshold changed: {key}")
        thresholds={fraction:float(np.quantile(fitted_presence,fraction)) for fraction in THRESHOLDS}
        train=np.vstack((data_p,data_b))
        _,analog0=grid_validity_and_analog(scaled0,relevant0,train)
        _,analog21=grid_validity_and_analog(scaled21,relevant21,train)
        score0=np.full(GRID_N,np.nan,dtype=np.float32)
        score21=np.full(GRID_N,np.nan,dtype=np.float32)
        score0[relevant0_idx]=model_probs(clf,poly,scaled0[relevant0_idx])
        score21[relevant21_idx]=model_probs(clf,poly,scaled21[relevant21_idx])
        record.update(status="COMPLETE_05DEG_POTENTIAL_CLIMATIC_HABITAT",
                      model_reproduction_max_abs_error=err,
                      threshold_0BP_native10=threshold,
                      modern=count_available(free0>=5,finite0,analog0),
                      LGM=count_available(free21>=5,finite21,analog21))
        scenario_defs=[("primary_10pct_majority_4neighbor",.10,5,4),
                       ("neighbor8_10pct_majority",.10,5,8),
                       ("plant_05pct_majority",.05,5,4),
                       ("plant_15pct_majority",.15,5,4),
                       ("terrain_any_10pct",.10,1,4),
                       ("terrain_all_10pct",.10,9,4)]
        for name,q,terrain,near in scenario_defs:
            m0=suitable_mask(score0,thresholds[q],free0,analog0,finite0,
                             land_cut=terrain).reshape(GRID_SHAPE)
            m21=suitable_mask(score21,thresholds[q],free21,analog21,finite21,
                             land_cut=terrain).reshape(GRID_SHAPE)
            g0=graph_metrics(m0,lat,connectivity=near)
            g21=graph_metrics(m21,lat,connectivity=near)
            record[name]={"0BP":g0,"21ka":g21,
               "component_delta_21ka_minus_0BP":g21["components"]-g0["components"],
               "area_delta_km2":g21["suitable_km2"]-g0["suitable_km2"]}
            if name=="primary_10pct_majority_4neighbor":
                main_masks0.append((key,np.packbits(m0.reshape(-1))))
                main_masks21.append((key,np.packbits(m21.reshape(-1))))
        result.append(record);status[record["status"]]+=1
        if (j+1)%50==0:
            print(json.dumps({"source_plant_rows_processed":j+1,
                  "completed":status["COMPLETE_05DEG_POTENTIAL_CLIMATIC_HABITAT"],
                  "reproduction_error_max":max_source_error}),flush=True)
    if status["COMPLETE_05DEG_POTENTIAL_CLIMATIC_HABITAT"]!=416:
        raise ValueError("416 original source-qualified plant niche models failed to reproduce")
    if len(main_masks0)!=len(main_masks21)!=416:
        raise ValueError("main plant surface status count mismatch")
    maskids=[a[0] for a in main_masks0]
    if maskids!=[b[0] for b in main_masks21]:raise ValueError("0BP/21ka masks plant ordering drift")
    with (args.out/"all_518_plant_potential_2d_graph_v02.csv").open("w",encoding="utf8",newline="") as f:
        keys=sorted(set().union(*(r.keys() for r in result)))
        writer=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
        writer.writeheader()
        for x in result:
            writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,dict) else v
                             for k,v in x.items()})
    np.savez_compressed(args.out/"plant_potential_climate_graph_masks_primary_v02.npz",
                        plant_wcvp_ids=np.asarray(maskids,dtype="<U16"),
                        mask_0BP=np.stack([x[1] for x in main_masks0]),
                        mask_21ka=np.stack([x[1] for x in main_masks21]),
                        grid_shape=np.asarray(GRID_SHAPE,dtype=int),
                        lat_centers=np.asarray(lat,dtype=np.float32),
                        land_cut=np.array(5),threshold_quantile=np.array(.10),
                        nonanalog_excluded=np.array(True))
    modeled=[x for x in result if x["status"]=="COMPLETE_05DEG_POTENTIAL_CLIMATIC_HABITAT"]
    main=[x["primary_10pct_majority_4neighbor"] for x in modeled]
    d=np.array([v["area_delta_km2"] for v in main],float)
    source={"schema":"chocho_native_WCVP_416_plant_05deg_paleo_land_connectivity_v02",
        "status":"PHYSICAL_LAND_AND_CLIMATIC_POTENTIAL_GRAPH_NOT_LGM_OCCUPANCY_OR_GENETICS",
        "source_original_plant_models":416,
        "models_reproduced_on_original_6337_sites":416,
        "original_model_score_replay_max_abs_error":max_source_error,
        "original_presence_threshold_replay_max_abs_error":max_threshold_error,
        "source_grid_shape":list(GRID_SHAPE),"source_grid_cells":GRID_N,
        "plant_graphs_completed":len(modeled),
        "median_paired_potential_suitable_area_change_km2":float(np.median(d)),
        "plants_potential_suitable_area_decrease":int(np.sum(d<0)),
        "plants_potential_suitable_area_increase":int(np.sum(d>0)),
        "plants_potential_suitable_area_unchanged":int(np.sum(d==0)),
        "plant_graph_status_counts":dict(status),"main_source":{"land":">=5 of 9 PMIP4 10-minute icefree subcells",
         "climate":"CHELSA 0BP and 21ka, four original bioclimatic variables",
         "nonanalog":"Exclude modern training-envelope nonanalog cells at both epochs",
         "connectivity":"4-neighbor, dateline-connected, 0.5 degree",
         "threshold":"0BP plant native observed fit 10th percentile"},
        "source_sha256":dict(ROOT_SHA),
        "limitations":["Potential continuously adjacent climate-suitable land in a coarse 0.5° model; not occupied 21ka habitat",
            "Even moderate 0.5° PMIP4 land mask can create or erase narrow bridges and coastlines",
            "Model is trained on currently native documented WCVP regions and modern GBIF sampling",
            "Original 416 plant fits are replayed; 0.5° pixels are model predictions and are not independent plant observations",
            "Excluding climate-nonanalog pixels is a conservative uncertainty censoring convention, not proof of physical barriers",
            "No inferred biological migration, glacial plant occupancy, observed insect genetics or confirmatory host-specific genetic effect"]}
    (args.out/"summary.json").write_text(json.dumps(source,indent=2,sort_keys=True)+"\n")
    print(json.dumps(source,sort_keys=True))

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Species-specific native-source 0BP climate SDMs projected onto CHELSA 21ka.

Predicts relative *climate suitability* on 6343 source representative sites.
NOT actual LGM occupancy or an ecological/genetic dispersal corridor.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import roc_auc_score

PLAN_POINTS_SHA="f580dfc76c7eb3d2b7ffa68e54e35bf6534578465d40e997f5bcfec0f19efc79"
PLAN_CENSUS_SHA="1b009db9cc813a5ab992dd1baace23ea1aeed3673dbae5064a642be876c489c0"
SUPPORT_SHA="d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f"
NATIVE_UNITS_SHA="8fbf5e0fab0df919862955f725a1822c7af05e3d59887d785cbc3d22b9b96dc0"
TIME={"0BP":20,"LGM":-190}
VARS=("bio01","bio07","bio12","bio15")
THRESHOLDS=(.05,.10,.15)
MAX_BACKGROUND=1000

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def read_csv(p):
    with Path(p).open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def background_support_indices(taxon,native,support):
    by_region=defaultdict(list)
    for i,r in enumerate(support):by_region[r["area_code_l3"]].append(i)
    lookup=sorted({i for region in native.get(taxon,set()) for i in by_region.get(region,[])})
    if len(lookup)>MAX_BACKGROUND:
        lookup=sorted(lookup,key=lambda i:hashlib.sha256(
            f"chocho-native-background-v01|{taxon}|{i}".encode()).hexdigest())[:MAX_BACKGROUND]
    return np.array(sorted(lookup),dtype=int)

def climate_sample(paths,coords):
    import rasterio
    data=[]
    for p in paths:
        with rasterio.open(p) as src:
            if not src.crs or not src.crs.is_geographic:
                raise ValueError("CHELSA climate must be geographic")
            vals=np.array([float(x[0]) for x in src.sample(coords,masked=False)],float)
            if src.nodata is not None:vals[vals==float(src.nodata)]=np.nan
            data.append(vals)
    return np.column_stack(data)

def block_label(lat,lon):
    return (int(math.floor((float(lat)+90)/2)),int(math.floor((float(lon)+180)/2)))

def spatial_holdout_labels(plant_coords,bg_coords,taxon):
    blocks=[block_label(p[0],p[1]) for p in np.vstack([plant_coords,bg_coords])]
    unique=sorted(set(blocks))
    if len(unique)<5:return np.zeros(len(blocks),bool)
    keyed=sorted(unique,key=lambda b:hashlib.sha256(
        f"chocho-native-spatial-fold-v01|{taxon}|{b[0]}|{b[1]}".encode()).hexdigest())
    hold=set(keyed[::5])
    return np.array([block in hold for block in blocks],bool)

def fit_predict(pres,bg,pres_coord,bg_coord,global_current,global_lgm,taxon):
    # All four variables have already been standardized using source-global 0BP.
    pp=np.asarray(pres,float);bb=np.asarray(bg,float)
    if pp.ndim!=2 or bb.ndim!=2 or pp.shape[1]!=4 or bb.shape[1]!=4 or len(pp)<30 or len(bb)<30:
        raise ValueError("inadequate native plant sources or background")
    if not np.isfinite(pp).all() or not np.isfinite(bb).all():
        raise ValueError("nonfinite model predictor source")
    poly=PolynomialFeatures(degree=2,include_bias=False)
    x=np.vstack([pp,bb]);y=np.r_[np.ones(len(pp)),np.zeros(len(bb))]
    feats=poly.fit_transform(x)
    if feats.shape[1]!=14:raise ValueError("four-variable degree2 feature count drift")
    fit=LogisticRegression(C=.3,solver="liblinear",class_weight="balanced",
                           max_iter=1000,random_state=0,tol=1e-5)
    fit.fit(feats,y)
    hold=spatial_holdout_labels(pres_coord,bg_coord,taxon)
    cv_auc=None;cv_pres_n=int(sum(hold[:len(pp)]));cv_bg_n=int(sum(hold[len(pp):]))
    if (cv_pres_n>=5 and cv_bg_n>=5 and
        sum(~hold[:len(pp)])>=20 and sum(~hold[len(pp):])>=20):
        loo=LogisticRegression(C=.3,solver="liblinear",class_weight="balanced",
                               max_iter=1000,random_state=0,tol=1e-5)
        loo.fit(feats[~hold],y[~hold])
        cv_auc=float(roc_auc_score(y[hold],loo.predict_proba(feats[hold])[:,1]))
    score_pres=fit.predict_proba(poly.transform(pp))[:,1]
    now=fit.predict_proba(poly.transform(global_current))[:,1]
    then=fit.predict_proba(poly.transform(global_lgm))[:,1]
    lo=np.min(x,axis=0);hi=np.max(x,axis=0)
    margin=.05*np.maximum(hi-lo,.001)
    nonanalog=np.any((global_lgm<lo-margin)|(global_lgm>hi+margin),axis=1)
    threshold=np.quantile(score_pres,.10)
    current_suitable=now>=threshold
    lgm_suitable=then>=threshold
    analog=~nonanalog
    if not np.any(analog):raise ValueError("no analog 21ka support sites")
    return {
        "source_presence_n":len(pp),"source_background_n":len(bb),
        "cv_auc_presence_vs_background":cv_auc,
        "heldout_presence_n":cv_pres_n,"heldout_background_n":cv_bg_n,
        "threshold_0BP_presence_quantile_10":float(threshold),
        "climate_suitable_fraction_now":float(np.mean(current_suitable)),
        "climate_suitable_fraction_21ka":float(np.mean(lgm_suitable)),
        "climate_suitable_fraction_delta_21ka_minus_now":float(np.mean(lgm_suitable)-np.mean(current_suitable)),
        "nonanalog_21ka_site_fraction":float(np.mean(nonanalog)),
        "supported_analog_21ka_site_n":int(sum(analog)),
        "analog_only_suitable_fraction_now":float(np.mean(current_suitable[analog])),
        "analog_only_suitable_fraction_21ka":float(np.mean(lgm_suitable[analog])),
        "sensitivity_05_10_15":{str(q):{
            "present":float(np.mean(now>=np.quantile(score_pres,q))),
            "LGM":float(np.mean(then>=np.quantile(score_pres,q)))} for q in THRESHOLDS}
    },now,then,nonanalog

def main():
    args=argparse.ArgumentParser()
    for key in ("plant-plan","thinned-native","native-units","support","chelsa-root","out"):
        args.add_argument("--"+key,type=Path,required=True)
    a=args.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for p,digest in ((a.plant_plan,PLAN_CENSUS_SHA),(a.thinned_native,PLAN_POINTS_SHA),
                     (a.native_units,NATIVE_UNITS_SHA),(a.support,SUPPORT_SHA)):
        if sha(p)!=digest:raise ValueError("source changed before LGM hindcast: "+str(p))
    plants=read_csv(a.plant_plan);support=read_csv(a.support)
    if len(plants)!=518 or len(support)!=6373:raise ValueError("WCVP source dimensions drift")
    regions=defaultdict(set)
    for r in read_csv(a.native_units):
        regions[str(int(r["accepted_plant_name_id"]))].add(r["area_code_l3"])
    pres=defaultdict(list)
    with gzip.open(a.thinned_native,"rt",encoding="utf8",newline="") as f:
        for r in csv.DictReader(f):pres[r["wcvp_id"]].append(r)
    if sum(map(len,pres.values()))!=73386:raise ValueError("source 10km-native-presence sample size drift")
    # All climate samples source by exact current presence locations and the 6373 WCVP support positions.
    coords=[];lookup={}
    def ix(rows):
        result=[]
        for lon,lat in rows:
            p=(round(float(lon),8),round(float(lat),8))
            if p not in lookup:lookup[p]=len(coords);coords.append(p)
            result.append(lookup[p])
        return np.array(result,dtype=int)
    global_idx=ix([(r["longitude"],r["latitude"]) for r in support])
    source_idx={}
    for r in plants:
        key=str(int(r["wcvp_id"]))
        if r["status"]=="READY_NATIVE_SDM":
            source_idx[key]=ix([(p["longitude"],p["latitude"]) for p in pres[key]])
    if len(source_idx)!=416:raise ValueError("416 WCVP-native eligible plant taxa did not reconcile")
    current_paths=[a.chelsa_root/f"CHELSA_TraCE21K_{v}_20_V1.0.tif" for v in VARS]
    lgm_paths=[a.chelsa_root/f"CHELSA_TraCE21K_{v}_-190_V1.0.tif" for v in VARS]
    for p in current_paths+lgm_paths:
        if not p.is_file():raise FileNotFoundError(str(p))
    print(json.dumps({"native_plant_climate_point_locations":len(coords),
                      "global_support_locations":len(global_idx)}),flush=True)
    current=climate_sample(current_paths,coords)
    past_raw=climate_sample(lgm_paths,[coords[i] for i in global_idx])
    source_ref=current[global_idx]
    valid_now=np.isfinite(source_ref).all(axis=1)
    if sum(valid_now)!=6343:raise ValueError("original 6343 global present climate support drift")
    mean=np.mean(source_ref[valid_now],axis=0)
    sd=np.std(source_ref[valid_now],axis=0)
    if np.any(sd<=0):raise ValueError("climate variable standard deviation is zero")
    current=(current-mean)/sd
    past=(past_raw-mean)/sd
    now_global=current[global_idx]
    comparable=valid_now & np.isfinite(past).all(axis=1)
    if sum(comparable)<6000:raise ValueError("too few source-fixed comparable global reference sites")
    global_now=now_global[comparable];global_lgm=past[comparable]
    ids=[str(int(r["wcvp_id"])) for r in plants]
    results=[];now_predictions=np.full((len(plants),sum(comparable)),np.nan,dtype=np.float32)
    lgm_predictions=np.full_like(now_predictions,np.nan)
    model_result=Counter()
    bg_table=defaultdict(list)
    for j,r in enumerate(support):bg_table[r["area_code_l3"]].append(j)
    for k,plant in enumerate(plants):
        taxon=ids[k];base={"wcvp_id":taxon,"accepted_host_name":plant["accepted_host_name"],
                            "native_regions":int(plant["native_region_L3_count"]),
                            "original_native_occurrence_n":int(plant["original_native_observations"]),
                            "native_thinned_n":int(plant["thinned_native_observations"]),
                            "background_n":int(plant["native_background_support_points"])}
        if plant["status"]!="READY_NATIVE_SDM":
            base["model_status"]=plant["status"];results.append(base);model_result[plant["status"]]+=1
            continue
        bg_idx=sorted({i for reg in regions[taxon] for i in bg_table[reg]})
        if len(bg_idx)>MAX_BACKGROUND:
            bg_idx=sorted(bg_idx,key=lambda i:hashlib.sha256(
                f"chocho-native-background-v01|{taxon}|{i}".encode()).hexdigest())[:MAX_BACKGROUND]
        bg_idx=np.array(sorted(bg_idx),dtype=int)
        if len(bg_idx)!=int(plant["native_background_support_points"]):
            raise ValueError("native region background selection drift")
        pi=source_idx[taxon]
        p_env=current[pi];b_env=now_global[bg_idx]
        pxy=np.array([(float(x["latitude"]),float(x["longitude"])) for x in pres[taxon]])
        bxy=np.array([(float(support[i]["latitude"]),float(support[i]["longitude"])) for i in bg_idx])
        fp=np.isfinite(p_env).all(axis=1);fb=np.isfinite(b_env).all(axis=1)
        if sum(fp)<30 or sum(fb)<30:
            base["model_status"]="NONFINITE_NATIVE_MODEL_TRAINING"
            base["finite_presence_n"]=int(sum(fp));base["finite_background_n"]=int(sum(fb))
            results.append(base);model_result[base["model_status"]]+=1
            continue
        try:
            fit,pred0,pred21,nonanalog=fit_predict(
                p_env[fp],b_env[fb],pxy[fp],bxy[fb],global_now,global_lgm,taxon)
            base.update(fit,model_status="COMPLETE_POTENTIAL_CLIMATE_HINDCAST")
            now_predictions[k]=pred0.astype(np.float32)
            lgm_predictions[k]=pred21.astype(np.float32)
        except (ValueError,OverflowError) as err:
            base.update(model_status="MODEL_NONIDENTIFIABLE",error=str(err)[:240])
        results.append(base);model_result[base["model_status"]]+=1
        if (k+1)%50==0:
            print(json.dumps({"plants_processed":k+1,"total":518,
                              "models_complete":model_result["COMPLETE_POTENTIAL_CLIMATE_HINDCAST"]}),flush=True)
    if len(results)!=518:raise RuntimeError("lost WCVP source species")
    keys=sorted(set().union(*(r.keys() for r in results)))
    with (a.out/"all_518_native_host_hindcast.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader()
        for r in results:
            row={k:(json.dumps(v,sort_keys=True) if isinstance(v,dict) else v) for k,v in r.items()}
            w.writerow(row)
    # Store native host-specific modern/LGM model scores across identical global support sites.
    np.savez_compressed(a.out/"native_host_climate_surface_scores.npz",
         plant_wcvp_id=np.array(ids),global_support_indices=np.where(comparable)[0],
         support_longitude=np.array([float(support[i]["longitude"]) for i in np.where(comparable)[0]],np.float32),
         support_latitude=np.array([float(support[i]["latitude"]) for i in np.where(comparable)[0]],np.float32),
         scores_present=now_predictions,scores_21ka=lgm_predictions)
    complete=[r for r in results if r["model_status"]=="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"]
    crossvalidated=[r["cv_auc_presence_vs_background"] for r in complete
                    if r["cv_auc_presence_vs_background"] is not None]
    summary={"schema":"chocho_WCVP_native_plant_TraCE21ka_potential_suitability_v0.1",
       "status":"EXPLORATORY_CLIMATIC_SUITABILITY_NOT_ACTUAL_LGM_PLANT_OCCUPANCY",
       "plant_source_taxa":518,"source_native_paired_insect_species":392,
       "input_model_eligible_taxa":416,"model_status_counts":dict(model_result),
       "model_complete_taxa":len(complete),
       "global_reference_sites_present_finite":int(sum(valid_now)),
       "global_sites_comparable_0BP_vs_21ka":int(sum(comparable)),
       "median_delta_lgm_minus_modern_global_site_fraction":float(np.median(
            [r["climate_suitable_fraction_delta_21ka_minus_now"] for r in complete])) if complete else None,
       "median_21ka_nonanalog_site_fraction":float(np.median(
            [r["nonanalog_21ka_site_fraction"] for r in complete])) if complete else None,
       "spatial_auc_available":len(crossvalidated),
       "median_spatial_auc_presence_vs_background":float(np.median(crossvalidated)) if crossvalidated else None,
       "climate_sources":["CHELSA-TraCE21k original 0BP and 21ka, bio01/07/12/15"],
       "no_lgm_occupancy_or_genetic_response_claim":True,
       "limits":["Presence-background logistic scores are not prevalence-calibrated occupancy probabilities",
           "Current WCVP L3 native regions define training background, not ancestral glacial regions",
           "Global model predictions are on 6343 support samples, not area-weighted full raster maps",
           "No dispersal, historical ice, sea-level/shoreline or observed gene flow modeled",
           "Plant GBIF iNaturalist observational effort is nonrandom",
           "Spatial holdout AUC is a screening diagnostic not a validated independent distribution prediction"]}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))

if __name__=="__main__":main()

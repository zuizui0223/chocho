#!/usr/bin/env python3
"""Prediction gate for 225 butterflies × 3,696 islands, NOT port causality.

Compare:
 M0: historical candidate-butterfly reporting, island area, latitude/longitude,
     and source butterfly species fixed effects;
 M1: M0 + *modern indexed WPI port* geography;
 M2: M1 + global introduced-host resource portfolio × WPI port association.

Five-fold SPATIAL-group holdout at 10°x10° groups, then 2018–2022 temporal
rank validation on still-unindexed species-island pairs. The GEB host portfolio
is a global resource opportunity, NOT historically present alien plants on the
specific island. Modern port locations are not 2001–2010 ship connectivity.
Do not interpret index probabilities as true occupancy, new arrival or breeding.
"""
from __future__ import annotations
import csv,gzip,hashlib,json,math
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss,average_precision_score,brier_score_loss,roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import OneHotEncoder,StandardScaler

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
PATH=HERE/"global_225x3696_butterfly_island_risk.csv.gz"
SRC=HERE/"global_original_14m_source_and_risk_receipt.json"
PORTS=HERE/"global_WPI_ports_assigned_to_GSHHG.csv"
OUT=HERE/"global_port_host_spatial_prediction_results.json"
COVARS=("log_area","log_pre_plant_exposure_proxy","pre_richness","abs_lat",
        "sin_lon","cos_lon","lat2","lon2")
# Here log_pre_plant_exposure_proxy actually means logged *butterfly*
# observations, not plant data. Legacy name is avoided in user-visible results.

def read_input():
    if not PATH.is_file() or not SRC.is_file():
        raise ValueError("Global public source-rebuilt butterfly risk file missing; run 14m original global GIS before modeling")
    receipt=json.loads(SRC.read_text(encoding="utf8"))
    if receipt["status"]!="ORIGINAL_PUBLIC_GBIF_DERIVED_GLOBAL_ISLAND_RISK_NOT_CAUSAL_BUTTERFLY_COLONIZATION":
        raise ValueError("Global original GBIF source results not successfully archived")
    if receipt["global_all_butterfly_species_x_island_risk_cells"]!=831600:
        raise ValueError("Full source 225×3696 denominator not verified")
    df=pd.read_csv(PATH,compression="gzip",low_memory=False,
                   dtype={"species":str,"speciesKey":str,"landmass_id":str})
    if len(df)!=831600:
        raise ValueError("Missing butterfly island risk cells")
    if df.groupby(["speciesKey","landmass_id"]).size().max()!=1:
        raise ValueError("Duplicate species-island risk identity")
    if not set(df["eligible_2011_2017"].dropna().unique())<={0,1}:
        raise ValueError("Invalid first period risk")
    if not set(df["first_2011_2017"].dropna().unique())<={0,1}:
        raise ValueError("Invalid first period endpoint")
    if any((df["first_2011_2017"]>df["eligible_2011_2017"])
            | (df["first_2018_2022"]>df["eligible_2018_2022"])):
        raise ValueError("Future indexed event outside first-risk cohort")
    return df,receipt

def transformed(df):
    x=pd.DataFrame(index=df.index)
    area=np.maximum(df["GSHHG_equal_area_km2"].to_numpy(dtype=float),.001)
    pre=np.maximum(df["pre2011_total_candidate_butterfly_GBIF_records"].to_numpy(dtype=float),0)
    richness=np.maximum(df["pre2011_island_candidate_species_index_count"].to_numpy(dtype=float),0)
    lat=df["reference_latitude"].to_numpy(dtype=float)
    lon=df["reference_longitude"].to_numpy(dtype=float)
    x["log_area"]=np.log1p(area)
    x["log_pre_plant_exposure_proxy"]=np.log1p(pre)
    x["pre_richness"]=np.log1p(richness)
    x["abs_lat"]=abs(lat)/90
    x["sin_lon"]=np.sin(np.deg2rad(lon))
    x["cos_lon"]=np.cos(np.deg2rad(lon))
    x["lat2"]=(lat/90)**2
    x["lon2"]=(lon/180)**2
    x["modern_WPI_any"]=df["source_2019_WPI_port_count"].gt(0).astype(float).to_numpy()
    x["log_modern_WPI_ports"]=np.log1p(np.maximum(df["source_2019_WPI_port_count"].to_numpy(dtype=float),0))
    host=df["GEB_introduced_share_of_contemporary"].to_numpy(dtype=float)
    missing=~np.isfinite(host)
    if missing.any():
        # Missing GEB match: never fill with zero biological opportunity.
        x["host_data_available"]=~missing
        x["host_fraction"]=np.where(missing,np.nanmedian(host),host)
    else:
        x["host_data_available"]=1
        x["host_fraction"]=host
    x["host_x_modern_WPI"]=x["modern_WPI_any"]*x["host_fraction"]
    x["host_x_log_port_count"]=x["log_modern_WPI_ports"]*x["host_fraction"]
    if any(not np.isfinite(x[c]).all() for c in x.columns):
        raise ValueError("Missing or nonfinite regression feature")
    x["species"]=df["speciesKey"].astype(str)
    # Landscape spatial blocks avoid hiding adjacent island replicates in
    # opposite crossvalidation folds (10-degree approximate, not archipelago).
    y=np.floor((lat+90)/10).astype(int)
    z=np.floor((lon+180)/10).astype(int)
    blocks=y*36+z
    return x,blocks,missing.sum()

def fit_compare(train,test,y,models):
    """Source-fold-only scaling + species indicators; report actual held-out scores."""
    enc=OneHotEncoder(handle_unknown="ignore",sparse_output=True)
    sp_train=enc.fit_transform(train[["species"]])
    sp_test=enc.transform(test[["species"]])
    outputs={}
    for name,cols in models.items():
        ss=StandardScaler()
        x_train=ss.fit_transform(train[cols])
        x_test=ss.transform(test[cols])
        a=sparse.hstack([sp_train,sparse.csr_matrix(x_train)],format="csr")
        b=sparse.hstack([sp_test,sparse.csr_matrix(x_test)],format="csr")
        model=LogisticRegression(C=1.0,max_iter=320,solver="lbfgs")
        model.fit(a,y[train.index])
        pred=np.clip(model.predict_proba(b)[:,1],1e-10,1-1e-10)
        yy=y[test.index]
        outputs[name]={"log_loss":float(log_loss(yy,pred,labels=[0,1])),
            "brier":float(brier_score_loss(yy,pred)),
            "average_precision":float(average_precision_score(yy,pred)) if yy.sum() else None,
            "positive_rate":float(yy.mean()),"positive_events":int(yy.sum()),
            "n":len(yy),"predictors":cols}
    return outputs

def run():
    raw,source=read_input()
    data,block,missing=transformed(raw)
    early=raw["eligible_2011_2017"].to_numpy(dtype=int)==1
    df=data.loc[early].reset_index(drop=True)
    blocks=block[early]
    y=raw.loc[early,"first_2011_2017"].to_numpy(dtype=int)
    if len(y)<600000 or y.sum()<50:
        raise ValueError("Unexpected global early risk support")
    names={
       "M0_historical_butterfly_index_area_geography_species":list(COVARS),
       "M1_plus_modern_WPI_ports":list(COVARS)+["modern_WPI_any","log_modern_WPI_ports"],
       "M2_plus_species_introduced_host_port_interaction":
            list(COVARS)+["modern_WPI_any","log_modern_WPI_ports",
                          "host_data_available","host_x_modern_WPI",
                          "host_x_log_port_count"]
    }
    groups=np.unique(blocks)
    if len(groups)<15:raise ValueError("Spatial crossvalidation blocks absent")
    folds=[]
    split=GroupKFold(n_splits=5)
    for k,(itr,itest) in enumerate(split.split(df,y,groups=blocks),1):
        # train/test groups cannot overlap, response cannot be used for choice
        # of islands/plant predictor.
        if set(blocks[itr])&set(blocks[itest]):
            raise ValueError("Geographically leaked holdout")
        train=df.iloc[itr].copy();test=df.iloc[itest].copy()
        train.index=np.arange(len(train))
        test.index=np.arange(len(test))+len(train)
        yy=np.concatenate([y[itr],y[itest]])
        fits=fit_compare(train,test,yy,names)
        folds.append({"fold":k,"n_training":len(itr),
          "n_test":len(itest),"n_spatial_training_blocks":len(set(blocks[itr])),
          "n_spatial_holdout_blocks":len(set(blocks[itest])),
          "results":fits})
        print("GLOBAL_FOLD",k,"positive",int(y[itest].sum()),
              "M0_M1_M2_BRIER",*[round(fits[m]["brier"],7) for m in names],flush=True)
    summary={}
    for model in names:
        fold_values=[f["results"][model] for f in folds]
        n=np.array([v["n"] for v in fold_values],float)
        summary[model]={
           "spatial_block_heldout_sample_weighted_logloss":float(np.average([v["log_loss"] for v in fold_values],weights=n)),
           "spatial_block_heldout_sample_weighted_brier":float(np.average([v["brier"] for v in fold_values],weights=n)),
           "spatial_block_heldout_unweighted_AP":float(np.mean([v["average_precision"] for v in fold_values if v["average_precision"] is not None])),
           "per_fold_brier":[v["brier"] for v in fold_values],
           "per_fold_index_events":[v["positive_events"] for v in fold_values]
        }
    b0=summary["M0_historical_butterfly_index_area_geography_species"]["spatial_block_heldout_sample_weighted_brier"]
    b1=summary["M1_plus_modern_WPI_ports"]["spatial_block_heldout_sample_weighted_brier"]
    b2=summary["M2_plus_species_introduced_host_port_interaction"]["spatial_block_heldout_sample_weighted_brier"]
    report={
      "schema":"chocho_global_butterfly_new_record_port_resource_spatial_CV_v1",
      "status":"GLOBAL_OBSERVATIONAL_PREDICTION_NOT_CAUSAL_HOST_ESTABLISHMENT",
      "checked_at_utc":datetime.now(timezone.utc).isoformat(),
      "original_global_GBIF_doi":source["original_doi"],
      "input_global_butterfly_island_225_risk_cells":len(raw),
      "early_2011_2017_eligible_cells":int(early.sum()),
      "first_indexed_2011_2017":int(y.sum()),
      "GEB_host_metric_missing_data_cells_across_full_panel":int(missing),
      "spatial_crossvalidation":"5 GroupKFold folds of island 10x10 latitude-longitude degree blocks",
      "source_block_count":len(groups),
      "models":summary,
      "full_source_fold_details":folds,
      "delta_brier_M1_vs_M0":b1-b0,
      "delta_brier_M2_vs_M1":b2-b1,
      "relative_brier_improvement_M1_vs_M0":(b0-b1)/b0 if b0 else None,
      "relative_brier_improvement_M2_vs_M1":(b1-b2)/b1 if b1 else None,
      "endpoint":"new indexed GBIF occurrence at previously unindexed island 2011-2017",
      "historical_detection_effort_proxy":"GBIF indexed pre2011 butterfly observations and richness; incomplete effort correction",
      "GEB_host_resource_proxy":"species global fraction of contemporary potential host resource in introduced regions; NOT actual island host index",
      "modern_2019_WPI_proxy":"cartographic harbor presence and counts, NOT 2001-2010 maritime cargo",
      "limitations":[
        "Not an island-arrival, original introduction, establishment, natural dispersal or ship-mediated causal model.",
        "No island-specific historically introduced host plants for all 3696 global islands.",
        "Unobserved butterfly index and unindexed plants do NOT imply ecological absence.",
        "Causal port effects are confounded by island area, climate, human density, recording effort and colonial trade histories even when held-out prediction improves.",
        "Island spatial folds are 10-degree geographic cells and are not independent geological archipelagos or ocean-basin holdouts.",
        "The species-wide GEB resource portfolio and butterfly records share taxonomy and historical sources.",
        "Original 225 source butterfly species do not represent every global Lepidoptera species.",
        "Modern 2019 WPI island port assignments include 130 ambiguous WPI coastal islands; sensitivity pending."
      ]
    }
    OUT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:report[k] for k in (
        "status","input_global_butterfly_island_225_risk_cells",
        "early_2011_2017_eligible_cells","first_indexed_2011_2017",
        "relative_brier_improvement_M1_vs_M0","relative_brier_improvement_M2_vs_M1")},
        ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":run()

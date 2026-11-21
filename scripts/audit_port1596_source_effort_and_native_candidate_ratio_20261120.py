#!/usr/bin/env python3
"""Exploratory country-stratified observation-effort diagnostic for 1596 WPI ports.

Source 1: original COMPLETE 1596-port cached HOSTS-catalogue status counts,
on original WPI port IDs (original Git blob fdf73b11...).
Source 2: independent NEW 2026-11-20 GBIF index, 1596 ports, datasetKey facets,
with source WPI centres, but NOT original historical WKT buffered data.

Both taxon-status classes use same original plant-GBIF sampling; this is not
field host quality/abundance, causal adjustment, or the independent WCVP126
port-versus-whole-introduced-flora main paper result. Do not retroactively
assign this new source to original GBIF historical providers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

ORIGINAL_PORTS = Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
NEW_PROVIDER_SOURCE = Path("results/geb_v05_gbif_source_attribution_20261120/NEW_GBIF_ORIGINAL1596_PORTS_DATASETKEY_INDEX_20261120.json")
ORIGINAL_BLOB = "fdf73b11ddcb1644f46730123350ae52cc61cae6"
NEW_SHA256 = "e79c3ad437d484424d41a5fc8e784c9c829070d740f1120be769d5a2d69548b6"

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def source(original: Path, current: Path) -> pd.DataFrame:
    raw = (original / ORIGINAL_PORTS).read_bytes()
    original_blob = hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if original_blob != ORIGINAL_BLOB:
        raise RuntimeError("Original 1596-port historical source identity drift")
    original_rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    old=pd.DataFrame(original_rows)
    if len(old)!=3630 or old["port_id"].nunique()!=3630:
        raise RuntimeError("Original original-WPI port roster drift")
    if sum(old.query_complete.eq("True"))!=1596:
        raise RuntimeError("Original 1596 source-eligible port cohort changed")
    data=(current/NEW_PROVIDER_SOURCE).read_bytes()
    if digest(data)!=NEW_SHA256:
        raise RuntimeError("New 1596 indexed GBIF provider aggregate is not source-exact")
    new=json.loads(data)
    if (new["schema"]!="chocho_all1596_original_WPI2019_ports_NEW_current_GBIF_datasetKey_facets_20261120_v2"
        or not new["audit_no_missing_duplicate_port_ID"] or
        not new["audit_individual_port_and_four_shard_provider_count_exact"] or
        new["source_original_2026_GBIF_contributing_datasetKeys_recovered"]):
        raise RuntimeError("Current GBIF index wrongly relabelled historical or source unverifiable")
    readings=[]
    for row in new["new_source_port_index_facets"]:
        providers=row["current_GBIF_record_counts_by_datasetKey"]
        count=int(row["new_2026_GBIF_index_count"])
        if sum(int(p["count"]) for p in providers)!=count:
            raise RuntimeError("Port source provider record sums do not match total")
        readings.append({"port_id":str(row["original_WPI2019_OBJECTID"]),
            "current_GBIF_records":count,
            "current_GBIF_datasets":int(row["new_datasetKey_count"]),
            "top1_current_source_share":max([int(p["count"]) for p in providers],default=0)/max(1,count)})
    modern=pd.DataFrame(readings)
    if len(modern)!=1596 or modern.port_id.nunique()!=1596:
        raise RuntimeError("New GBIF re-index source has missing/duplicate original port")
    sample=old[old.query_complete.eq("True")].merge(modern,how="left",on="port_id",validate="one_to_one")
    if len(sample)!=1596 or sample[["current_GBIF_records","current_GBIF_datasets"]].isna().any().any():
        raise RuntimeError("Original-to-new port source join lost rows")
    for col in ["contemporary_not_native_in_catalogue", "native_in_catalogue",
                "ambiguous_across_regions","unknown_catalogue_status"]:
        sample[col]=pd.to_numeric(sample[col],errors="raise")
        if (sample[col]<0).any():
            raise RuntimeError("Original WCVP/HOSTS status candidate counts invalid")
    sample["country"] = sample["country"].replace("",pd.NA).fillna("UNKNOWN")
    return sample

def estimate(df: pd.DataFrame, response: str, columns: list[str]) -> dict:
    design=df[["country","native_in_catalogue","contemporary_not_native_in_catalogue",
               "current_GBIF_records","current_GBIF_datasets"]].copy()
    design["candidate_log"]=np.log1p(design.contemporary_not_native_in_catalogue)
    design["native_log"]=np.log1p(design.native_in_catalogue)
    design["log_candidate_native_ratio"]=np.log((design.contemporary_not_native_in_catalogue+0.5)
                                                 /(design.native_in_catalogue+0.5))
    design["log_current_GBIF_records"]=np.log1p(design.current_GBIF_records)
    design["log_current_GBIF_datasets"]=np.log1p(design.current_GBIF_datasets)
    names=[response]+columns
    within=design[names]-design.groupby("country")[names].transform("mean")
    fit=sm.OLS(within[response],within[columns]).fit(
        cov_type="cluster",cov_kwds={"groups":design.country,"use_correction":True})
    if fit.params.isna().any() or not np.isfinite(fit.conf_int().to_numpy()).all():
        raise RuntimeError("Source country clustered regression not estimable")
    ci=fit.conf_int()
    return {"n_ports":len(df),"country_clusters":int(df.country.nunique()),
            "response":response,"country_fixed_effects":True,
            "cluster_robust_error":"country CR1 statsmodels",
            "within_R2":float(fit.rsquared),
            "regressors":{x:{"beta":float(fit.params[x]),"se":float(fit.bse[x]),
                              "p":float(fit.pvalues[x]),"CI95":[float(ci.loc[x,0]),float(ci.loc[x,1])]}
                          for x in columns}}

def audit(df: pd.DataFrame) -> dict:
    candidates="contemporary_not_native_in_catalogue"
    native="native_in_catalogue"
    records="current_GBIF_records"
    providers="current_GBIF_datasets"
    rho={}
    for out in [candidates,native]:
        for inp in [records,providers]:
            r=spearmanr(df[out],df[inp])
            rho[out+"__"+inp]={"rho":float(r.statistic),
                "p":float(r.pvalue),"p_is_float_underflow_zero":bool(r.pvalue==0)}
    source_predictor=spearmanr(df[records],df[providers])
    single=["log_current_GBIF_records"]
    both=["log_current_GBIF_records","log_current_GBIF_datasets"]
    primary={a:estimate(df,a,single) for a in ["candidate_log","native_log","log_candidate_native_ratio"]}
    joint={a:estimate(df,a,both) for a in ["candidate_log","native_log","log_candidate_native_ratio"]}
    subsets={
        "current_GBIF_positive_ports":df[df.current_GBIF_records>0],
        "both_original_catalogue_status_classes_positive":
            df[(df[candidates]>0)&(df[native]>0)],
        "original_country_groups_at_least_3_ports":
            df[df.country.map(df.country.value_counts())>=3],
        "exclude_current_one_provider_dominance_90pct":
            df[df.top1_current_source_share<0.9]
    }
    sensitivity={k:{"single_record_model":estimate(q,"log_candidate_native_ratio",single),
                    "joint_record_and_publisher_model":estimate(q,"log_candidate_native_ratio",both)}
                  for k,q in subsets.items()}
    loo=[]
    for group in sorted(df.country.unique()):
        s=df[df.country!=group]
        design=pd.DataFrame({"country":s.country.to_numpy(),
              "response":np.log((s[candidates].to_numpy()+0.5)/(s[native].to_numpy()+0.5)),
              "x":np.log1p(s.current_GBIF_records.to_numpy())})
        z=design[["x","response"]]-design.groupby("country")[["x","response"]].transform("mean")
        beta=(z["x"]*z["response"]).sum()/(z["x"]**2).sum()
        loo.append(float(beta))
    if len(loo)!=df.country.nunique() or not np.isfinite(loo).all():
        raise RuntimeError("Country leave-one-out slopes invalid")
    current_positive=df[df.current_GBIF_records>0]
    out={
      "schema":"chocho_GEBv05_source_observation_effort_candidate_native_port1596_20261120_v1",
      "status":"EXPLORATORY_CROSS_SNAPSHOT_SOURCE_EFFORT_DIAGNOSTIC_NOT_HISTORIC_PROVENANCE_OR_MAIN_WCVP_RESULT",
      "original_source_git_blob":ORIGINAL_BLOB,
      "new_GBIF_2026_source_sha256":NEW_SHA256,
      "n_orig_WPI_complete_ports":len(df),"country_groups_including_UNKNOWN":int(df.country.nunique()),
      "ports_missing_original_country":int(df.country.eq("UNKNOWN").sum()),
      "new_GBIF_current_positive_ports":len(current_positive),
      "new_GBIF_current_zero_index_ports":int(df.current_GBIF_records.eq(0).sum()),
      "original_candidate_not_native_catalogue_zero_ports":int(df[candidates].eq(0).sum()),
      "original_native_catalogue_zero_ports":int(df[native].eq(0).sum()),
      "new_current_provider_top1_share_median_positive_ports":
          float(current_positive.top1_current_source_share.median()),
      "new_current_provider_top1_share_more_than_0_5_all_ports":
          int(df.top1_current_source_share.gt(.5).sum()),
      "new_current_provider_top1_share_at_least_0_9_all_ports":
          int(df.top1_current_source_share.ge(.9).sum()),
      "cross_snapshot_Spearman_associations":rho,
      "new_GBIF_records_and_datasets_Spearman_rho":float(source_predictor.statistic),
      "country_cluster_primary_record_associations":primary,
      "country_cluster_two_predictor_sensitivity":joint,
      "country_cluster_log_candidate_native_ratio_sensitivities":sensitivity,
      "leave_one_original_country_out_log_ratio_record_gradient":{
        "N_dropped_countries":len(loo),"min_beta":min(loo),"max_beta":max(loo),
        "count_negative":sum(x<0 for x in loo),
        "is_not_equivalence_proof":True
      },
      "ecological_qualification":[
        "Original contemporary-not-native catalogue status is a candidate source classification, NOT independently documented establishment through ports.",
        "Absolute candidate and native status richness strongly correlates with public GBIF botanical recording in both source classes.",
        "Country-clustered ratio model fails to reject zero, but does not prove ecological equality or absence of effort bias.",
        "The current 20-Nov GBIF source count and provider count are highly collinear; any marginal multivariable provider coefficient is post hoc and sensitivity-dependent.",
        "New 2026 sampled source is from a different time and geoDistance area than original 2026 speciesKey WKT port cache; do not use these regressions to retrospectively validate original cached datasetKey publisher identities.",
        "Most crucially, the MAIN 126-botanical-island WCVP whole-flora versus WPI port-count regression has a different outcome/source panel and cannot be refuted or confirmed by this auxiliary 1596-port index diagnostic.",
        "The original 1596 sampling is partial (3630 WPI candidates) and geographically ordered; country fixed effects are insufficient to estimate true complete-port biodiversity or causal introduction.",
        "All country CIs are cluster robust descriptive in reused ports. This exploratory diagnostic was specified after examining raw rho; do not label preregistered inference."
      ]}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--current-root",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    x=p.parse_args()
    s=source(x.original_root,x.current_root)
    result=audit(s)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("WPI1596_PORT_INTRO_NATIVE_SOURCE_EFFORT_COUNTRY_FE_QA_PASS",
        json.dumps({"raw_rho_introduced":result["cross_snapshot_Spearman_associations"]["contemporary_not_native_in_catalogue__current_GBIF_records"]["rho"],
        "raw_rho_native":result["cross_snapshot_Spearman_associations"]["native_in_catalogue__current_GBIF_records"]["rho"],
        "ratio_beta":result["country_cluster_primary_record_associations"]["log_candidate_native_ratio"]["regressors"]["log_current_GBIF_records"]["beta"],
        "ratio_CI95":result["country_cluster_primary_record_associations"]["log_candidate_native_ratio"]["regressors"]["log_current_GBIF_records"]["CI95"]},ensure_ascii=False))
if __name__=="__main__":
    main()

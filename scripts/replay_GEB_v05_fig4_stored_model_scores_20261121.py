#!/usr/bin/env python3
"""Source-exact but limited Fig4 endpoint SCORE replay for original GBIF cohorts.

Recalculate species-equal 155-index forecast improvement from original
published model error table and region-equal 224-species assemblage
improvement from all original 169 held-out region errors. This does not
re-fit classifiers, recompute held-out predictions or remake original
GBIF event/model input from raw occurrence records. Post-hoc 29-draw
thinning/world transfer summaries require different original outputs
and are NOT independently reproduced here.

Both independent endpoints share original GBIF DOI 10.15468/dl.umqmsf,
which is not historical 225x51 or 1596-port query provenance.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path

SOURCE={
 "forecast_metrics":"3e2dcf109b48d3c3766e34ff7603d5061cb850b2",
 "region_holdouts":"be2018b452e6c44eb3c60282e9d4af6414d9a70c",
}
MODES={
 "temporal":(208482,4.120175607931824),
 "species":(208482,-0.2085201539930192),
 "spatial":(208482,2.9952059516988134),
 "spatial_buffer100":(102728,2.7356579110032984),
 "temporal_effort20":(159084,2.975749318968573),
 "temporal_radius250":(180654,3.7626496823350153),
 "temporal_exclude_top_provider":(191004,5.545564256849527),
 "temporal_precision100":(135532,5.344940058277037)
}
REGION_MODELS=("M0_all_flora_effort_prior","M1_known_host","M2_other_aliens","M3_both")
EXPECTED_REGIONAL_MSE={
 "M0_all_flora_effort_prior":0.0038319659427122565,
 "M1_known_host":0.0038259328683706347,
 "M2_other_aliens":0.003831744481720401,
 "M3_both":0.003828246452242267
}
def read(path,pin):
    data=path.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if blob!=pin:raise AssertionError("Original GBIF derived source score Git blob changed: "+path.name)
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
def replay(forecast,region):
    f=read(forecast,SOURCE["forecast_metrics"])
    if len(f)!=len(MODES)*7:raise AssertionError("Expected original 8 forecast modes x 7 models")
    forecasts={}
    for mode,(expected_rows,expected_gain) in MODES.items():
        r=[x for x in f if x["mode"]==mode]
        names={x["model"] for x in r}
        if len(r)!=7 or len(names)!=7:
            raise AssertionError("Forecast model group missing or duplicated "+mode)
        if any(int(x["rows"])!=expected_rows or int(x["species"])!=155 for x in r):
            raise AssertionError("Forecast source tested cohort drift "+mode)
        values={x["model"]:float(x["species_equal_mse"]) for x in r}
        reference=values["M3_both"]
        with_hosts=values["M5_host_composition"]
        if reference<=0 or with_hosts<=0:raise AssertionError("Invalid original forecast errors")
        gain=100*(reference-with_hosts)/reference
        if abs(gain-expected_gain)>1e-10:
            raise AssertionError("Source forecast percent gain drift "+mode)
        forecasts[mode]={
            "species":155,"evaluated_species_cell_cases":expected_rows,
            "reference_climate_landcover_M3_species_equal_MSE":reference,
            "plus_known_larval_host_resource_and_composition_M5_species_equal_MSE":with_hosts,
            "source_recomputed_MSE_gain_pct":gain,
            "model_refitted_from_raw_GBIF_occurrences":False
        }
    if set(x["mode"] for x in f)!=set(MODES):
        raise AssertionError("Unexpected original forecast modes")
    s=read(region,SOURCE["region_holdouts"])
    if len(s)!=4*169:
        raise AssertionError("Expected four 169-region held-out original butterfly score panels")
    scores=defaultdict(dict)
    level={}
    for row in s:
        model=row["model"]
        reg=row["wgsrpd3_code"]
        if model not in REGION_MODELS or reg in scores[model]:
            raise AssertionError("Unknown/duplicate 224-region original model and region score")
        if int(row["test_pairs"])!=168 or int(row["training_pairs"])!=14028:
            raise AssertionError("Dependent 169-region pairwise holdout population changed")
        mse=float(row["mse"])
        if not math.isfinite(mse) or mse<0:
            raise AssertionError("Unusable original 224-species region model held-out loss")
        scores[model][reg]=mse
        if reg in level and level[reg]!=row["world_division"]:
            raise AssertionError("Source region world-division changed")
        level[reg]=row["world_division"]
    if len(level)!=169 or any(set(t)!=set(level) for t in scores.values()):
        raise AssertionError("Four original region holdouts do not evaluate same 169 regions")
    means={k:sum(x.values())/len(x) for k,x in scores.items()}
    for k,target in EXPECTED_REGIONAL_MSE.items():
        if abs(means[k]-target)>1e-12:
            raise AssertionError("Original source region-equal MSE changed "+k)
    gain=100*(means["M2_other_aliens"]-means["M3_both"])/means["M2_other_aliens"]
    if abs(gain-0.09129078138747806)>1e-10:
        raise AssertionError("Original documented H-after-O improvement differs")
    folds_better=sum(scores["M3_both"][reg]<scores["M2_other_aliens"][reg] for reg in level)
    if folds_better!=96:
        raise AssertionError("Original 224-region H-after-O heldout regional win count changed")
    return {
      "schema":"GEB_v05_Fig4_two_GBIF_endpoints_score_table_replay_not_raw_model_fit_v1",
      "status":"TWO_ORIGINAL_155_AND_224_GBIF_ENDPOINT_LOSS_TABLES_RECALCULATED_NO_RETRAINING_OR_RECORDLEVEL_SOURCE",
      "GBIF_download_DOI":"10.15468/dl.umqmsf",
      "GBIF_original_224_butterfly_download_records":14637152,
      "two_model_endpoints_not_independent_data_downloads":True,
      "forecast155":forecasts,
      "observed_community224":{
        "n_species":224,"n_WGSRPD3_regions":169,
        "dependent_region_pairs":14196,"held_region_splits":169,
        "region_equal_MSE_by_model":means,
        "source_recomputed_known_host_H_after_other_aliens_O_MSE_gain_pct":gain,
        "regions_where_M3_both_better_than_M2_other_aliens":folds_better,
        "source_input_is_original_169_region_model_holdout_MSE_not_individual_GBIF_occurrences":True
      },
      "other_Fig4_source_claims_not_recomputed":[
        "29-draw original datasetKey-day/record-total matched thinning responses",
        "eight-world held-region transfer re-fits and South-American dyad calibration",
        "host plant species source extraction, predictor age filtering and model training",
        "butterfly occurrence/absence/detection or causal host-resource effect"
      ],
      "source_scope":"Original GBIF model output errors re-aggregated only, not fully independent observed community/regression replay",
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--forecast-metrics",type=Path,required=True)
    p.add_argument("--observed224-region-losses",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    report=replay(a.forecast_metrics,a.observed224_region_losses)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GEB_FIG4_155_8_FORECAST_MODES_SCORE_SOURCE_RECALCULATED_PASS")
    print("GEB_FIG4_224_169_HELD_REGION_H_AFTER_O_SCORE_SOURCE_RECALCULATED_PASS")
    print("GEB_FIG4_RAW_GBIF_MODEL_RETRAINING_AND_THINNING_NOT_CLAIMED")
    print(json.dumps({"temporal_155_gain":report["forecast155"]["temporal"]["source_recomputed_MSE_gain_pct"],
       "species_unseen_155_gain":report["forecast155"]["species"]["source_recomputed_MSE_gain_pct"],
       "spatial_155_gain":report["forecast155"]["spatial"]["source_recomputed_MSE_gain_pct"],
       "distance100_155_gain":report["forecast155"]["spatial_buffer100"]["source_recomputed_MSE_gain_pct"],
       "community224_H_after_O_gain":report["observed_community224"]["source_recomputed_known_host_H_after_other_aliens_O_MSE_gain_pct"]},ensure_ascii=False))
if __name__=="__main__":
    main()

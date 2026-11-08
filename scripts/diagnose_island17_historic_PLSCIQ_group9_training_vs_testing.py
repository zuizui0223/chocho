#!/usr/bin/env python3
"""POST-RESULT source-only decomposition of Borneo original group9 influence.

No new fitted models: use only the previously archived original 2048 historical
PLSCIQ two-island-unseen source-pair predictions.
"""
from __future__ import annotations
import argparse,csv,gzip,json,math
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/"results/island17_historic_PLSCIQ_LO2_20261108"
SOURCE=PARENT/"LO2_heldout_predictions.csv.gz"
PANEL=PARENT/"original17_136_source_pair_manifest.csv"
PARENT_SCORE=PARENT/"comparison_scope_summary.csv"
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
MODELS=("BASE","MED_SUM","MED_SUM_DIFF","MAX_SUM_DIFF")
COMPARE=(
    ("BASE_to_MED_SUM","BASE","MED_SUM"),
    ("BASE_to_MED_SUM_DIFF","BASE","MED_SUM_DIFF"),
    ("MED_SUM_to_MED_SUM_DIFF","MED_SUM","MED_SUM_DIFF"),
    ("BASE_to_MAX_SUM_DIFF","BASE","MAX_SUM_DIFF"),
    ("MED_SUM_DIFF_to_MAX_SUM_DIFF","MED_SUM_DIFF","MAX_SUM_DIFF"),
)
SCORES=(
    ("ALL17_ALL136","all17","all",136),
    ("FULL17_TRAINED_NO9_TEST120","all17","no9",120),
    ("REFIT_NO9_TEST120","leave_original_landmass9_out","no9",120),
    ("FULL17_TRAINED_GROUP9_TEST16","all17","touch9",16)
)

def read(path):
    with Path(path).open("rt",newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def read_gzip(path):
    with gzip.open(path,"rt",newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f)

def write(path,rows):
    with Path(path).open("wt",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def near(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10)

def run(outdir):
    manifest=read(PANEL)
    if len(manifest)!=136 or len({int(p["original_pair_source_index"]) for p in manifest})!=136:
        raise ValueError("Original complete source 17-island/136-pair sample missing")
    include9={int(x["original_pair_source_index"]) for x in manifest
              if x["original_island_a"]=="9" or x["original_island_b"]=="9"}
    omit9={int(x["original_pair_source_index"]) for x in manifest
           if x["original_island_a"]!="9" and x["original_island_b"]!="9"}
    if len(include9)!=16 or len(omit9)!=120 or bool(omit9&include9):
        raise ValueError("Source group9 16/120 pair partition changed")
    original_scores=read(PARENT_SCORE)
    if len(original_scores)!=20:
        raise ValueError("Original frozen PLSCIQ score table missing")
    lookup={}
    fullsets=defaultdict(set)
    for row in read_gzip(SOURCE):
        scope,target,model=row["scope"],row["target"],row["model"]
        index=int(row["original_pair_source_index"])
        key=(scope,target,model,index)
        if key in lookup or target not in TARGETS or model not in MODELS:
            raise ValueError("Source double-unseen prediction duplicate/missing model")
        if (scope=="all17" and index not in (omit9|include9) or
            scope=="leave_original_landmass9_out" and index not in omit9):
            raise ValueError("Unrecognized saved original source target/pair population")
        if scope not in ("all17","leave_original_landmass9_out"):
            raise ValueError("Unrecognized source model training population")
        y=float(row["observed"]);predict=float(row["predicted"])
        if not 0<=y<=1 or not 0<=predict<=1:
            raise ValueError("Historical source original Jaccard or model forecasts outside 0..1")
        sq=(y-predict)**2;ae=abs(y-predict)
        if not near(sq,row["squared_error"]) or not near(ae,row["absolute_error"]):
            raise ValueError("Saved original historical source model squared or absolute error changed")
        lookup[key]=(y,predict,sq,ae)
        fullsets[(scope,target,model)].add(index)
    if len(lookup)!=2048:
        raise ValueError("All 2048 source original LO2 observations must remain")
    for scope in ("all17","leave_original_landmass9_out"):
        expected=(omit9|include9) if scope=="all17" else omit9
        for target in TARGETS:
            for model in MODELS:
                if fullsets[(scope,target,model)]!=expected:
                    raise ValueError("Missing original source model/scope island-pair identity")
    # Direct same original test pair comparison between the frozen 17 fit and 16 fit.
    for target in TARGETS:
        for model in MODELS:
            for index in omit9:
                parent=lookup[("all17",target,model,index)]
                child=lookup[("leave_original_landmass9_out",target,model,index)]
                if not near(parent[0],child[0]):
                    raise ValueError("Different actual butterfly outcome on same target/source pair")
    scores=[]
    for title,scope,where,n in SCORES:
        indices=sorted((omit9|include9) if where=="all" else omit9 if where=="no9" else include9)
        if len(indices)!=n:
            raise ValueError("Original source-only control/test stratum size differs")
        for target in TARGETS:
            for label,ref,aug in COMPARE:
                x=[lookup[(scope,target,ref,i)] for i in indices]
                y=[lookup[(scope,target,aug,i)] for i in indices]
                a=math.fsum(z[2] for z in x)/n
                b=math.fsum(z[2] for z in y)/n
                c=math.fsum(z[3] for z in x)/n
                d=math.fsum(z[3] for z in y)/n
                scores.append({
                    "evaluation_scope":title,"source_model_training_scope":scope,
                    "test_source_scope":where,"target":target,"contrast":label,
                    "original_source_unique_test_pairs":n,
                    "original_source_train_pairs_per_fold":105 if scope=="all17" else 91,
                    "MSE_reference":a,"MSE_augmented":b,
                    "MSE_reduction_percent":100*(a-b)/a,
                    "MAE_reference":c,"MAE_augmented":d,
                    "MAE_reduction_percent":100*(c-d)/c,
                    "pairs_lower_squared_error":sum(t[2]>v[2] for t,v in zip(x,y)),
                    "pairs_higher_squared_error":sum(t[2]<v[2] for t,v in zip(x,y)),
                    "only_saved_source_predictions_no_new_models":True,
                    "post_result_diagnostic_not_preregistered":True
                })
    if len(scores)!=40:
        raise ValueError("Expected exact 4 target strata x 2 butterfly periods x 5 contrasts")
    # Require prior 20 stored primary+no9 scores not changed.
    bysource={(r["scope"],r["target"],r["contrast"]):r for r in original_scores}
    for r in scores:
        if r["evaluation_scope"] not in ("ALL17_ALL136","REFIT_NO9_TEST120"):continue
        oldscope="all17" if r["evaluation_scope"]=="ALL17_ALL136" else "leave_original_landmass9_out"
        old=bysource[(oldscope,r["target"],r["contrast"])]
        for k in ("MSE_reference","MSE_augmented","MSE_reduction_percent",
                  "MAE_reference","MAE_augmented","MAE_reduction_percent"):
            if not near(r[k],old[k]):
                raise ValueError("Original 20 source-frozen historical PLSCIQ MSE/MAE comparisons changed: "+k)
    out=Path(outdir).resolve()
    if out.is_relative_to(ROOT/"results") or out.is_relative_to(ROOT/"provenance") or out.is_relative_to(ROOT/"data"):
        raise ValueError("Do not overwrite original source predictions")
    out.mkdir(parents=True,exist_ok=True)
    write(out/"group9_training_vs_testing_loss_decomposition_40.csv",scores)
    result={
        "schema":"chocho_historical_PLSCIQ_group9_original_training_testing_decomposition_postresult_v1",
        "status":"EXACT_ORIGINAL_2048_SAVED_FORECASTS_40_GROUP9_DECOMPOSITION_CONTRASTS_AUDITED",
        "original_source_full_islands":17,
        "original_source_no9_islands":16,
        "original_source_full_test_pairs":136,
        "source_pairs_not_containing_group9":120,
        "source_pairs_touching_group9":16,
        "source_predictions_reused_without_refit":len(lookup),
        "original_precommitted_20_comparisons_reproduced_exactly":True,
        "new_postresult_summary_rows":len(scores),
        "explanation":"score distribution vs model refit: compare FULL17_TRAINED_NO9_TEST120 with REFIT_NO9_TEST120",
        "no_independent_pair_level_pvalues_or_causal_mechanism":True
    }
    (out/"group9_diagnostic_receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("GROUP9_TRAIN_TEST_DIAGNOSTIC_JSON="+json.dumps(scores,sort_keys=True),flush=True)
    print("GROUP9_TRAIN_TEST_DIAGNOSTIC_COUNTS="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)

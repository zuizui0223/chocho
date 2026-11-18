#!/usr/bin/env python3
"""Independent standard-library replay of original group9 post-result 40 MSE/MAE scores.

Uses archived original prediction records and original source pair manifest only.
No import of diagnostic/fitting functions. Requires exact original 20 results.
"""
from __future__ import annotations
import argparse,csv,gzip,math,json
from collections import defaultdict
from pathlib import Path
PARENT=Path(__file__).resolve().parent
CASES={
 "ALL17_ALL136":("all17",None,136),
 "FULL17_TRAINED_NO9_TEST120":("all17",False,120),
 "REFIT_NO9_TEST120":("leave_original_landmass9_out",False,120),
 "FULL17_TRAINED_GROUP9_TEST16":("all17",True,16)
}
CONTRASTS={
 "BASE_to_MED_SUM":("BASE","MED_SUM"),
 "BASE_to_MED_SUM_DIFF":("BASE","MED_SUM_DIFF"),
 "MED_SUM_to_MED_SUM_DIFF":("MED_SUM","MED_SUM_DIFF"),
 "BASE_to_MAX_SUM_DIFF":("BASE","MAX_SUM_DIFF"),
 "MED_SUM_DIFF_to_MAX_SUM_DIFF":("MED_SUM_DIFF","MAX_SUM_DIFF")
}
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")

def read(f):
    with Path(f).open("r",encoding="utf-8",newline="") as s:
        return list(csv.DictReader(s))

def gz(f):
    with gzip.open(f,"rt",encoding="utf-8",newline="") as s:
        yield from csv.DictReader(s)

def check(a,b,key):
    if not math.isclose(float(a),float(b),abs_tol=1e-10,rel_tol=1e-10):
        raise ValueError("Independent group9 source decomposition numeric mismatch "+key)

def audit(out):
    manifest=read(PARENT/"original17_136_source_pair_manifest.csv")
    allid={int(r["original_pair_source_index"]):r for r in manifest}
    if len(allid)!=136:
        raise ValueError("Original historical index source pair manifest changed")
    table=defaultdict(dict)
    for r in gz(PARENT/"LO2_heldout_predictions.csv.gz"):
        group=(r["scope"],r["target"],r["model"])
        idx=int(r["original_pair_source_index"])
        if idx in table[group]:
            raise ValueError("Duplicate old model forecast")
        table[group][idx]=(float(r["squared_error"]),float(r["absolute_error"]))
    if sum(len(x) for x in table.values())!=2048:
        raise ValueError("Original source all 2048 double-island-unseen model predictions missing")
    published=read(Path(out)/"group9_training_vs_testing_loss_decomposition_40.csv")
    receipt=json.loads((Path(out)/"group9_diagnostic_receipt.json").read_text())
    if len(published)!=40 or receipt["original_precommitted_20_comparisons_reproduced_exactly"] is not True:
        raise ValueError("Original 40 source-only result rows/audit missing")
    seen=set()
    for r in published:
        case=r["evaluation_scope"];target=r["target"];contrast=r["contrast"]
        if case not in CASES or target not in TARGETS or contrast not in CONTRASTS:
            raise ValueError("Unexpected post-hoc source group9 comparison")
        ident=(case,target,contrast)
        if ident in seen:raise ValueError("Duplicate source group9 loss contrast")
        seen.add(ident)
        scope,touch,n=CASES[case]
        subset=sorted(i for i,item in allid.items()
             if touch is None or ((item["original_island_a"]=="9" or item["original_island_b"]=="9")==touch))
        if len(subset)!=n or int(r["original_source_unique_test_pairs"])!=n:
            raise ValueError("Source original group9 16 vs 120 pair test grouping changed")
        ref,aug=CONTRASTS[contrast]
        l0=[table[(scope,target,ref)][i] for i in subset]
        l1=[table[(scope,target,aug)][i] for i in subset]
        mse0=math.fsum(e[0] for e in l0)/n;mse1=math.fsum(e[0] for e in l1)/n
        mae0=math.fsum(e[1] for e in l0)/n;mae1=math.fsum(e[1] for e in l1)/n
        for key,value in (
            ("MSE_reference",mse0),("MSE_augmented",mse1),
            ("MSE_reduction_percent",100*(mse0-mse1)/mse0),
            ("MAE_reference",mae0),("MAE_augmented",mae1),
            ("MAE_reduction_percent",100*(mae0-mae1)/mae0),
            ("pairs_lower_squared_error",sum(x[0]>y[0] for x,y in zip(l0,l1))),
            ("pairs_higher_squared_error",sum(x[0]<y[0] for x,y in zip(l0,l1))),
        ):
            check(r[key],value,key)
        if r["only_saved_source_predictions_no_new_models"].lower()!="true" or r["post_result_diagnostic_not_preregistered"].lower()!="true":
            raise ValueError("Posthoc source observation was mislabeled as prespecified/causal")
    if len(seen)!=40 or receipt["new_postresult_summary_rows"]!=40:
        raise ValueError("Incomplete source group9 post-result evaluation panel")
    result={"schema":"chocho_independent_PLSCIQ_source_group9_train_test_decomposition_replay_v1",
        "status":"ORIGINAL_2048_LO2_FORECASTS_40_GROUP9_STRATA_SCORES_INDEPENDENTLY_RECONSTRUCTED",
        "original_source_model_predictions":2048,
        "source_original_pair_index_count":136,
        "original_group9_test_pairs":16,
        "original_no9_test_pairs":120,
        "source_decomposed_contrasts":40,
        "previously_frozen_models_not_refitted_or_selected":True,
        "posthoc_diagnostic_not_an_independent_causal_experiment":True}
    (Path(out)/"group9_independent_diagnostic_audit.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("GROUP9_INDEPENDENT_DIAGNOSTIC_AUDIT="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    audit(p.parse_args().out_dir)

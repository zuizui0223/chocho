#!/usr/bin/env python3
"""Independent stdlib audit of historical PLSCIQ original17 butterfly LO2 source.

Reconstructs original source WPI24 -> qualified22 -> 17 island source median/max,
original 136 source pair Jaccards, 136/120 both-endpoint-excluded LO2 folds,
all 2048 reported model predictions/losses and 20 MSE/MAE comparisons.
Does not import or trust the numerical fitting module; offline CI refits OLS.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,itertools,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
F17=ROOT/"results/island51_PortWatch17_geoverification_20261108/WPI17_independent_PortWatch_source_geography_LOCODE.csv"
F13=ROOT/"results/island51_PortWatch13_geoverification_20261108/WPI13_independent_PortWatch_source_geography_LOCODE.csv"
UNCTAD=ROOT/"results/island51_maritime_raw_recovery_20261108/UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv"
SOURCE=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
RECEIPT=ROOT/"results/island51_maritime_raw_recovery_20261108/UNCTAD_PLSCIQ_source_recovery_receipt.json"
WPI_SHA="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
MODELS=("BASE","MED_SUM","MED_SUM_DIFF","MAX_SUM_DIFF")
CONTRASTS={
 "BASE_to_MED_SUM":("BASE","MED_SUM"),
 "BASE_to_MED_SUM_DIFF":("BASE","MED_SUM_DIFF"),
 "MED_SUM_to_MED_SUM_DIFF":("MED_SUM","MED_SUM_DIFF"),
 "BASE_to_MAX_SUM_DIFF":("BASE","MAX_SUM_DIFF"),
 "MED_SUM_DIFF_to_MAX_SUM_DIFF":("MED_SUM_DIFF","MAX_SUM_DIFF"),
}
SCOPE=("all17","leave_original_landmass9_out")

def read(path):
    with Path(path).open("r",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def gzread(path):
    with gzip.open(path,"rt",newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f)

def near(a,b,label):
    if not math.isfinite(float(a)) or not math.isfinite(float(b)) or not math.isclose(
            float(a),float(b),rel_tol=2e-10,abs_tol=1e-10):
        raise ValueError(f"Historical PLSCIQ independent audit mismatch {label}: {a} vs {b}")

def audit(directory):
    out=Path(directory)
    raw=WPI.read_bytes()
    if hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()!=WPI_SHA:
        raise ValueError("Original WPI source geometry changed")
    wpi=read(WPI)
    if len(wpi)!=194 or len({r["island_id"] for r in wpi})!=51:
        raise ValueError("Original complete WPI194/51 source geometry missing")
    p17=read(F17);p13=read(F13);all_source=p17+p13
    if len(p17)!=17 or len(p13)!=13 or len({x["original_WPI_port_id"] for x in all_source})!=30:
        raise ValueError("Source original 30 candidate ports not preserved")
    verified={r["original_WPI_port_id"]:r for r in all_source if
        r["original_WPI_independent_geographic_identity_verified"]=="True"}
    if len(verified)!=24 or len({r["original_island_id"] for r in verified.values()})!=17:
        raise ValueError("24 source independent geographical WPI harbor IDs missing")
    indeces=read(UNCTAD)
    byoriginal={r["original_WPI_id"]:r for r in indeces}
    if len(indeces)!=194 or len(byoriginal)!=194:
        raise ValueError("Original historical UNCTAD source WPI194 page dataset incomplete")
    old=json.loads(RECEIPT.read_text())
    if (old["series"]!=1300 or old["total_port_quarter_rows"]!=81900 or
        old["nonmissing_port_quarter_rows"]!=57312 or
        not old["source_version"]["is_exact_2026_10_30_original_recovered"]):
        raise ValueError("The original UNCTAD 2026-10-30 edition cannot be independently verified")
    qualified={}
    for pid,r in verified.items():
        source=byoriginal[pid]
        if r["original_island_id"]!=source["original_island_id"]:
            raise ValueError("Wrong island for independent historical connectedness port")
        if (source["status"]=="UNIQUE_COUNTRY_AND_CANONICAL_NAME_CANDIDATE_NOT_GEO_VERIFIED" and
            int(source["available_2006_2010_quarters"] or 0)>=16):
            n=float(source["mean_2006_2010_PLCSIQ_old_index_base"])
            if not math.isfinite(n) or n<0:raise ValueError("Historical UNCTAD source index missing/nonfinite")
            qualified[pid]=(source["original_island_id"],n)
    if len(qualified)!=22:
        raise ValueError("Original historical 16/20-quarter source 22 WPI port support changed")
    med=defaultdict(list)
    for island,value in qualified.values():
        med[island].append(value)
    if len(med)!=17:
        raise ValueError("17 original historical connectedness island grouping changed")
    original_pairs=read(SOURCE)
    if len(original_pairs)!=1275:
        raise ValueError("Original complete 1275 butterfly pairs missing")
    orig={i:r for i,r in enumerate(original_pairs)}
    eligible={i for i,r in orig.items() if r["island_a"] in med and r["island_b"] in med}
    no9={i for i in eligible if orig[i]["island_a"]!="9" and orig[i]["island_b"]!="9"}
    if len(eligible)!=136 or len(no9)!=120:
        raise ValueError("Original source target-independent 17/16 island population is not C17,2/C16,2")
    savedports=read(out/"original24_geoverified_ports_historical_PLSCIQ_quality.csv")
    if len(savedports)!=24:
        raise ValueError("Original source port qualification ledger missing")
    for row in savedports:
        pid=row["original_WPI_port_id"]
        if pid not in verified or row["original_island_id"]!=verified[pid]["original_island_id"]:
            raise ValueError("Historical WPI index port ledger contains unverified source")
        src=byoriginal[pid]
        if row["available_2006_2010_quarters"]!=src["available_2006_2010_quarters"]:
            raise ValueError("Historical index port missingness changed")
        if (row["historical_index_eligible_geoverified_WPI_source_port"].lower()=="true")!=(pid in qualified):
            raise ValueError("Undated provider port promoted to historical data-ready")
        if src["mean_2006_2010_PLCSIQ_old_index_base"]:
            near(row["mean_2006_2010_PLCSIQ_old_index_base"],
                 src["mean_2006_2010_PLCSIQ_old_index_base"],"original old PLSCIQ mean")
    islands=read(out/"original51_historical_PLSCIQ_island_eligibility.csv")
    if len(islands)!=51 or {r["original_island_id"] for r in islands}!={r["island_id"] for r in wpi}:
        raise ValueError("Original all51 island eligibility and missingness ledger changed")
    for row in islands:
        island=row["original_island_id"]
        values=med.get(island,[])
        if ((row["historical_index_eligible"].lower()=="true")!=bool(values) or
            int(row["qualified_geoverified_historical_index_ports"])!=len(values)):
            raise ValueError("Source historical index island eligibility not observed")
        if values:
            near(row["historical_2006_2010_median_index"],statistics.median(values),"median source port index")
            near(row["historical_2006_2010_max_index"],max(values),"max source port index")
        elif row["historical_2006_2010_median_index"] or row["historical_2006_2010_max_index"]:
            raise ValueError("Missing indexed island falsely assigned zero historical index")
    pairframe=read(out/"original17_136_source_pair_manifest.csv")
    if len(pairframe)!=136:
        raise ValueError("Original 136 target-independent source pairs not preserved")
    seen=set()
    for row in pairframe:
        i=int(row["original_pair_source_index"])
        if i not in eligible or i in seen:
            raise ValueError("Source original butterfly pair index duplicated or outside historical island cohort")
        seen.add(i)
        p=orig[i]
        for f,sourcefield in (("original_island_a","island_a"),("original_island_b","island_b")):
            if row[f]!=p[sourcefield]:
                raise ValueError("Original butterfly historical island pair identity changed")
        for f in TARGETS:
            near(row[f],p[f],"observed original butterfly Jaccard")
        near(row["original_source_distance_km"],p["port_distance_km"],"original island source distance")
        if (row["sensitivity_original_group9_excluded"].lower()=="true")!=(i in no9):
            raise ValueError("Group9 source-only population freeze changed")
    if seen!=eligible:
        raise ValueError("Historical PLSCIQ original 136 eligible pairs missing")
    foldframe=read(out/"double_unseen_source_folds.csv")
    seen_folds=set()
    for row in foldframe:
        scope=row["scope"];i=int(row["original_pair_source_index"])
        dataset=eligible if scope=="all17" else no9 if scope=="leave_original_landmass9_out" else None
        if dataset is None or i not in dataset or (scope,i) in seen_folds:
            raise ValueError("Unknown/duplicated two-island-unseen original test pair")
        seen_folds.add((scope,i))
        p=orig[i]
        a,b=p["island_a"],p["island_b"]
        train={j for j in dataset if orig[j]["island_a"] not in (a,b) and orig[j]["island_b"] not in (a,b)}
        expected=105 if scope=="all17" else 91
        if (len(train)!=expected or int(row["train_pairs"])!=expected or
            int(row["test_pairs"])!=1 or row["test_island_a"]!=a or row["test_island_b"]!=b or
            row["both_test_endpoints_excluded_from_training"].lower()!="true"):
            raise ValueError("LO2 original test pair source endpoints included in train data")
    if len(foldframe)!=256 or len(seen_folds)!=256:
        raise ValueError("Original 136+120 pair folds incomplete")

    errors={}
    n=0
    for row in gzread(out/"LO2_heldout_predictions.csv.gz"):
        scope,target,model=row["scope"],row["target"],row["model"]
        i=int(row["original_pair_source_index"])
        key=(scope,target,model,i)
        if (scope not in SCOPE or target not in TARGETS or model not in MODELS or
            (scope,i) not in seen_folds or key in errors):
            raise ValueError("Original pair forecast unknown/duplicated")
        p=orig[i];actual=float(p[target]);estimate=float(row["predicted"])
        if not 0<=estimate<=1:
            raise ValueError("Predicted original Jaccard outside 0-1")
        near(row["observed"],actual,"original source true outcome")
        near(row["squared_error"],(actual-estimate)**2,"squared source model forecast error")
        near(row["absolute_error"],abs(actual-estimate),"absolute source model forecast error")
        if (row["test_island_a"]!=p["island_a"] or row["test_island_b"]!=p["island_b"] or
            int(row["train_pairs"])!=(105 if scope=="all17" else 91) or int(row["test_pairs"])!=1):
            raise ValueError("Original two-island-unseen forecast scope was changed")
        errors[key]=((actual-estimate)**2,abs(actual-estimate))
        n+=1
    if len(errors)!=2048 or n!=2048:
        raise ValueError("Original 2048 heldout forecasts not present exactly once")
    sums=read(out/"comparison_scope_summary.csv")
    if len(sums)!=20:
        raise ValueError("Expected exact 20 previously frozen source model contrasts")
    checks=set()
    for row in sums:
        scope,target,comparison=row["scope"],row["target"],row["contrast"]
        key=(scope,target,comparison)
        if scope not in SCOPE or target not in TARGETS or comparison not in CONTRASTS or key in checks:
            raise ValueError("Source 17-island model contrast unexpected/duplicated")
        checks.add(key)
        ref,aug=CONTRASTS[comparison]
        indices=sorted(eligible if scope=="all17" else no9)
        referr=[errors[(scope,target,ref,i)] for i in indices]
        augerr=[errors[(scope,target,aug,i)] for i in indices]
        rse=math.fsum(z[0] for z in referr)/len(referr)
        ase=math.fsum(z[0] for z in augerr)/len(augerr)
        rae=math.fsum(z[1] for z in referr)/len(referr)
        aae=math.fsum(z[1] for z in augerr)/len(augerr)
        if rse<=0 or rae<=0:raise ValueError("Historical source baseline has zero loss denominator")
        for field,value in (
            ("MSE_reference",rse),("MSE_augmented",ase),
            ("MSE_reduction_percent",100*(rse-ase)/rse),
            ("MAE_reference",rae),("MAE_augmented",aae),
            ("MAE_reduction_percent",100*(rae-aae)/rae),
            ("pairs_with_lower_squared_error",sum(a[0]>b[0] for a,b in zip(referr,augerr))),
            ("pairs_with_higher_squared_error",sum(a[0]<b[0] for a,b in zip(referr,augerr))),
            ("pairs_with_lower_absolute_error",sum(a[1]>b[1] for a,b in zip(referr,augerr))),
        ):
            near(row[field],value,"historical connectedness source contrast "+field)
    if len(checks)!=20:
        raise ValueError("Incomplete original historical source contrast set")
    conclusion={
      "schema":"chocho_independent_old_PLSCIQ_17source_island_butterfly_LO2_source_and_score_audit_v1",
      "status":"OLD_UNCTAD_SOURCE_22_PORTS_17_ISLANDS_2048_LO2_PREDICTIONS_ALL_VERIFIED",
      "original_source_WPI_ports":194,
      "original_51_islands_preserved":len(islands),
      "source_geographically_verified_index_ports":len(verified),
      "old_index_qualified_port_count":len(qualified),
      "unique_original_eligible_islands":len(med),
      "original_pair_target_outcomes_check":len(pairframe),
      "original_group9_excluded_pairs_check":len(no9),
      "original_double_unseen_fold_count":len(seen_folds),
      "all_individual_source_holdout_predictions_audited":n,
      "all_model_contrast_target_scope_results_audited":len(checks),
      "no_pair_level_p_value_inference":True,
      "historical_index_not_observed_ship_voyage_or_causal_colonization":True,
    }
    (out/"independent_source_and_prediction_audit.json").write_text(
        json.dumps(conclusion,sort_keys=True,indent=2)+"\n")
    print("HISTORICAL_PORT_PLSCIQ_INDEPENDENT_AUDIT="+json.dumps(conclusion,sort_keys=True),flush=True)
    return conclusion

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",required=True,type=Path)
    audit(p.parse_args().results_dir)

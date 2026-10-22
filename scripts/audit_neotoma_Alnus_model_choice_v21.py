#!/usr/bin/env python3
"""Exact source-age-model choice sensitivity for 21 source-labelled Alnus pollen samples.

No new LGM occurrence selection. Keep every original candidate, query no live
data, never treat default Neotoma chronology label as independently calibrated
or the genus-labelled pollen as observed occupied Alnus glutinosa trees.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED={
 "candidates":"aae646cfc66b7efc299a21b4a5ae83f52a9be90543ede2e5f50bedab63b5dbc8",
 "models":"4f957787237867d210b97bebda8a2ef581b8188ed21146d3419241719e43f386",
 "v20_summary":"6f69c5189746ab3fa4bfd7a2bde929828f15c15cb34aba4e78b05bc2b206c0df"
}
CAL="calibrated radiocarbon years bp"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def readcsv(path):
    with Path(path).open("r",encoding="utf8",newline="") as file:
        return list(csv.DictReader(file))

def numeric(value):
    if value is None or isinstance(value,bool) or str(value).strip()=="":
        return None
    try:x=float(value)
    except (ValueError,TypeError,OverflowError):return None
    return x if math.isfinite(x) else None

def yes(value):
    if value in (True,"True","true",1):return True
    if value in (False,"False","false",0):return False
    raise ValueError("source boolean field not recognized")

def choose_calibrated_default(row,models):
    original=str(row["source_chronologyid"])
    if original not in models:raise ValueError("source chronology not in 11 original models")
    source_model=models[original]
    if source_model["source_datasetid"]!=str(row["datasetid"]):
        raise ValueError("original pollen dataset does not match assigned chronology")
    source_cal=str(row["source_age_type"]).casefold()==CAL
    if not source_cal:raise ValueError("original v18 source labelled candidate no longer calibrated")
    age=numeric(row["sample_center_age_BP"])
    if age is None:raise ValueError("missing candidate center sample BP age")
    if yes(row["original_chronology_isdefault"])!=yes(source_model["isdefault"]):
        raise ValueError("source v20 original default chronology status changed")
    if yes(row["original_chronology_isdefault"]):
        default_model_id=original
        default_age=age
        default_type=row["source_age_type"]
        default_alternate=False
    else:
        alternates=json.loads(row["alternative_original_source_age_models"])
        # Keep only a default calibrated model on the EXACT source sampleid.
        valid=[x for x in alternates
           if yes(x.get("default_model_flag_from_chronology",False))
           and str(x.get("original_age_type","")).casefold()==CAL
           and numeric(x.get("original_age_numeric")) is not None]
        if len(valid)!=1:
            raise ValueError("source candidate lacks exactly one explicit calibrated default age model")
        default_model_id=str(valid[0]["chronologyid"])
        default_age=numeric(valid[0]["original_age_numeric"])
        default_type=valid[0]["original_age_type"]
        default_alternate=True
        if default_model_id not in models or not yes(models[default_model_id]["isdefault"]):
            raise ValueError("alternative source chronology falsely claimed default")
    m=models[default_model_id]
    if m["source_datasetid"]!=str(row["datasetid"]):
        raise ValueError("default chronology data source is not source pollen dataset")
    young=numeric(m["source_reliable_younger_BP"])
    old=numeric(m["source_reliable_older_BP"])
    if young is not None and old is not None and young>old:raise ValueError("invalid chronology source reliable bounds")
    inside=young<=default_age<=old if old is not None and young is not None else None
    earliest=numeric(m["source_control_depth_min_cm"])
    last=numeric(m["source_control_depth_max_cm"])
    sample_depth=numeric(row["sample_depth_cm"])
    in_depth=earliest<=sample_depth<=last if earliest is not None and last is not None and sample_depth is not None else None
    return {
        "datasetid":row["datasetid"],"sampleid":row["sampleid"],"siteid":row["siteid"],
        "positive_Alnus_pollen_datum_value":row["original_Alnus_pollen_positive_datum_value"],
        "original_selected_chronologyid":original,
        "original_selected_chronology_was_default":yes(row["original_chronology_isdefault"]),
        "original_selected_source_age_BP":age,
        "default_calibrated_chronologyid":default_model_id,
        "default_chronology_alternate_for_this_original_sample":default_alternate,
        "default_calibrated_source_age_BP":default_age,
        "default_calibrated_age_type":default_type,
        "default_calibrated_center_inside_18_24ka":18000<=default_age<=24000,
        "default_chronology_reliable_younger_BP":young,
        "default_chronology_reliable_older_BP":old,
        "default_center_within_reported_reliable_span":inside,
        "sample_depth_cm":sample_depth,
        "default_control_depth_min_cm":earliest,
        "default_control_depth_max_cm":last,
        "sample_depth_in_default_returned_control_extent":in_depth,
        "original_individual_sample_has_age_uncertainty_bounds":yes(row["source_sample_has_both_age_bounds"]),
        "default_model_control_count":int(m["source_age_controls_n"]),
        "metadata_only_candidate_not_independently_calibrated":True,
        "occupied_21ka_Alnus_population_confirmed":False
    }

def evaluate(candidate_rows,model_rows,source_v20):
    models={r["chronologyid"]:r for r in model_rows}
    if (len(candidate_rows)!=21 or len(models)!=11 or
       len(set((r["datasetid"],r["sampleid"]) for r in candidate_rows))!=21
       or source_v20["original_v18_source_labelled_18_24ka_positive_sample_n"]!=21):
        raise ValueError("21 unique source samples or full 11 source chronologies altered")
    output=[choose_calibrated_default(r,models) for r in candidate_rows]
    count=Counter()
    for x in output:
        count["source_original_labelled_18_24ka_samples"]+=1
        count["source_original_nondefault_to_default_calibrated_changed"]+=x["default_chronology_alternate_for_this_original_sample"]
        count["default_calibrated_center_inside_18_24ka"]+=x["default_calibrated_center_inside_18_24ka"]
        count["default_calibrated_center_outside_18_24ka"]+=not x["default_calibrated_center_inside_18_24ka"]
        count["default_calibrated_center_in_18_24ka_and_reported_reliable_span"]+=(
            x["default_calibrated_center_inside_18_24ka"] and x["default_center_within_reported_reliable_span"] is True)
        count["default_calibrated_center_18_24ka_both_reliable_and_controls_extent"]+=(
            x["default_calibrated_center_inside_18_24ka"] and x["default_center_within_reported_reliable_span"] is True
            and x["sample_depth_in_default_returned_control_extent"] is True)
        count["source_missing_individual_uncertainty_bounds"]+=not x["original_individual_sample_has_age_uncertainty_bounds"]
    by_dataset=defaultdict(list)
    for x in output:by_dataset[x["datasetid"]].append(x)
    site=[{"datasetid":did,"original_labelled_candidates":len(data),
           "default_calibrated_age_within_18_24ka":sum(x["default_calibrated_center_inside_18_24ka"] for x in data),
           "source_quality_metadata_triage_retained":sum(
             x["default_calibrated_center_inside_18_24ka"] and
             x["default_center_within_reported_reliable_span"] is True and
             x["sample_depth_in_default_returned_control_extent"] is True for x in data)}
          for did,data in sorted(by_dataset.items(),key=lambda x:int(x[0]))]
    summary={
     "schema":"chocho_Alnus_21ka_same_source_sample_default_chronology_choice_v21",
     "status":"EXPLORATORY_CHRONOLOGY_CHOICE_TRIAGE_NOT_CONFIRMED_21KA_LOCAL_ALNUS_OCCUPANCY",
     "sample_counts":dict(count),"by_source_pollen_dataset":site,
     "source_full_sample_n":21,"original_source_chronology_models_n":11,
     "reported_reliable_and_depth_control_metadata_only":True,
     "accepted_independently_calibrated_21ka_samples":0,
     "accepted_species_identified_local_Alnus_glutinosa_trees":0,
     "genetic_response_read":False,
     "source_sha256":EXPECTED,
     "limitations":["Age model default selection is a source metadata comparison, not an independent posterior probability for the true fossil age",
       "Source reliable age span and control depths do not replace missing sample-specific calibrated uncertainties",
       "Source positive Alnus glutinosa label might reflect transported Alnus pollen rather than local living host trees",
       "Original TTF/chocho recorded-host-versus-matched-other-real-plant historical specificity remains negative and unchanged"]
    }
    return output,summary

def main():
    a=argparse.ArgumentParser()
    for p in ("v20-samples","v20-models","v20-summary","out"):
        a.add_argument("--"+p,required=True,type=Path)
    z=a.parse_args();z.out.mkdir(parents=True,exist_ok=True)
    for key,path in (("candidates",z.v20_samples),("models",z.v20_models),("v20_summary",z.v20_summary)):
        if sha(path)!=EXPECTED[key]:raise ValueError("source v20 11 chronology models/21 samples SHA changed")
    rows,q=evaluate(readcsv(z.v20_samples),readcsv(z.v20_models),
                    json.loads(z.v20_summary.read_text()))
    path=z.out/"all_21_Alnus_original_vs_default_ages_v21.csv"
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    (z.out/"summary.json").write_text(json.dumps(q,indent=2,sort_keys=True)+"\n")
    print(json.dumps(q,sort_keys=True))

if __name__=="__main__":main()

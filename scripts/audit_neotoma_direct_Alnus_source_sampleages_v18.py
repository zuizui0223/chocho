#!/usr/bin/env python3
"""Direct Neotoma source sample pollen-datum and age-model inventory.

Read six SHA-pinned OFFICIAL full dataset downloads. Taxonid1299 Alnus glutinosa
datum is only a database classification; even with positive count and source
calibrated-age label it is not yet species-diagnostic occupied LGM vegetation.
No use of the unreliable occurrenceid -> sampleid endpoint as a join.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict

TAXON_ID=1299
CAL_TYPE="calibrated radiocarbon years bp"
STATUS="SOURCE_DATUM_SAMPLE_CHRONOLOGY_CANDIDATES_NOT_VERIFIED_LGM_OCCUPANCY"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def source_samples(raw,expected_dataset,expected_site):
    root=json.loads(raw)
    data=root.get("data") if isinstance(root,dict) else root
    if not isinstance(data,list) or len(data)!=1 or not isinstance(data[0],dict):
        raise ValueError("original full data download must have one source site")
    site=data[0].get("site")
    if not isinstance(site,dict) or str(site.get("siteid"))!=str(expected_site):
        raise ValueError("actual full data Neotoma site ID not expected frozen site")
    unit=site.get("collectionunit")
    dataset=unit.get("dataset") if isinstance(unit,dict) else None
    if not isinstance(dataset,dict) or str(dataset.get("datasetid"))!=str(expected_dataset):
        raise ValueError("full source data dataset ID not expected source")
    if dataset.get("datasettype") not in (None,"pollen"):
        raise ValueError("not source pollen data")
    samples=dataset.get("samples")
    if not isinstance(samples,list) or not all(isinstance(s,dict) for s in samples):
        raise ValueError("original full download lacks typed sample list")
    return samples

def finite_num(value):
    if isinstance(value,bool) or value is None:return None
    try:
        v=float(value)
    except (ValueError,TypeError,OverflowError):
        return None
    return v if math.isfinite(v) else None

def source_pollen_datum(sample):
    datums=sample.get("datum")
    if not isinstance(datums,list):raise ValueError("missing Neotoma sample datum array")
    matched=[v for v in datums if isinstance(v,dict) and str(v.get("taxonid"))==str(TAXON_ID)]
    positives=[v for v in matched if finite_num(v.get("value")) is not None and
               finite_num(v["value"])>0]
    zeros=[v for v in matched if finite_num(v.get("value"))==0]
    bad=[v for v in matched if finite_num(v.get("value")) is None or finite_num(v["value"])<0]
    return {"source_taxon_1299_datum_entries":len(matched),
            "source_positive_numeric_datums":len(positives),
            "source_zero_datums":len(zeros),
            "source_nonnumeric_or_negative_datums":len(bad),
            "total_positive_numeric_datum_value":sum(finite_num(v["value"]) for v in positives)}

def classify_source_age(age,contract):
    """A LABELLED candidate is NOT independently calibrated/verified chronology."""
    if not isinstance(age,dict):raise ValueError("invalid source age record")
    target_age=finite_num(age.get("age"))
    type_value=age.get("agetype")
    type_str=str(type_value or "")
    type_norm=" ".join(type_str.lower().split())
    chron=age.get("chronologyid")
    chron=str(int(chron)) if isinstance(chron,int) and not isinstance(chron,bool) else (
        chron.strip() if isinstance(chron,str) and chron.strip().isdigit() else None)
    lo=finite_num(age.get("ageyounger"))
    hi=finite_num(age.get("ageolder"))
    known_bounds=lo is not None and hi is not None and lo<=hi
    if lo is not None and hi is not None and lo>hi:
        raise ValueError("original source sample age uncertainty has inverted bounds")
    calibrated_label=type_norm==CAL_TYPE
    in_age=target_age is not None and contract["source"]["age_min_calibrated_BP"]<=target_age<=contract["source"]["age_max_calibrated_BP"]
    label_candidate=bool(calibrated_label and in_age and chron is not None)
    return {"source_age_BP_numeric":target_age,"source_agetype_label":type_str,
        "source_chronologyid":chron,"source_ageolder_BP":hi,"source_ageyounger_BP":lo,
        "has_both_source_age_bounds":known_bounds,
        "calibrated_age_label_not_independently_verified":calibrated_label,
        "source_age_within_18_24ka":bool(in_age),
        "source_labelled_18_24ka_cal_BP_candidate":label_candidate,
        "age_record_status":(
          "LABELLED_CAL_18_24KA_SOURCE_CANDIDATE_ONLY" if label_candidate else
          "OTHER_CHRONOLOGY_OR_TIME_SOURCE_ONLY")}

def validate_source(contract,data):
    byid=contract["source"]["dataset_sha256"]
    site=contract["source"]["dataset_site"]
    if (contract["status"]!="PREDECLARED_SAMPLE_LEVEL_SOURCE_CHECK_NO_ACCEPTED_LGM_OCCUPANCY"
        or set(byid)!=set(site) or len(byid)!=6
        or contract["source"]["taxonid"]!=TAXON_ID):
        raise ValueError("source fixed dataset and taxon identity contract drift")
    totals=Counter();source_samples_rows=[];source_ages_rows=[];summaries=[]
    seen_dataset_sample_ids=set()
    for did in sorted(byid,key=int):
        gz=data/f"dataset_{did}.json.gz"
        with gzip.open(gz,"rb") as fh:raw=fh.read()
        if sha(raw)!=byid[did]:raise ValueError(f"Neotoma v16 source dataset {did} SHA changed")
        samples=source_samples(raw,did,site[did])
        s_counts=Counter();used_ages=0;chronology_ids=set()
        for ordinal,sample in enumerate(samples):
            pollen=source_pollen_datum(sample)
            if not pollen["source_taxon_1299_datum_entries"]:continue
            sid=sample.get("sampleid")
            if sid is None:raise ValueError("direct source Alnus datum missing sampleid; do not invent ID")
            sid=str(sid)
            if (did,sid) in seen_dataset_sample_ids:
                raise ValueError("source dataset contains duplicate sampleid")
            seen_dataset_sample_ids.add((did,sid))
            s_counts["taxon1299_any_datum_sample_n"]+=1
            ages=sample.get("ages")
            if not isinstance(ages,list):raise ValueError("source sample ages absent")
            status="SOURCE_POSITIVE_DATUM_UNVERIFIED_POLLEN_SPECIES" if pollen["source_positive_numeric_datums"] else "ZERO_OR_UNKNOWN_DATUM_NOT_PRESENCE"
            if pollen["source_positive_numeric_datums"]:s_counts["positive_datum_sample_n"]+=1
            if len(ages)>1:s_counts["taxon1299_samples_multiple_age_models_n"]+=1
            row={"datasetid":did,"siteid":site[did],"sampleid":sid,
                 "analysisunitid":sample.get("analysisunitid"),
                 "sample_depth":sample.get("depth"),
                 "original_age_rows_n":len(ages),
                 "datum_status":status,
                 "source_species_label_not_morphology_verified":contract["source"]["source_taxon_label"],
                 "calibrated_18_24ka_sample_accepted_as_LGM_occupancy":False,
                 **pollen}
            candidate_count=0
            for age_i,age in enumerate(ages):
                item=classify_source_age(age,contract)
                if item["source_chronologyid"]:chronology_ids.add(item["source_chronologyid"])
                used_ages+=1
                if item["source_labelled_18_24ka_cal_BP_candidate"] and pollen["source_positive_numeric_datums"]:
                    candidate_count+=1
                    s_counts["labelled_candidate_age_model_rows_n"]+=1
                    if not item["has_both_source_age_bounds"]:
                        s_counts["candidate_missing_both_age_bounds_n"]+=1
                source_ages_rows.append({
                    "datasetid":did,"siteid":site[did],"sampleid":sid,
                    "source_age_record_ordinal":age_i,"datum_positive":bool(pollen["source_positive_numeric_datums"]),
                    "source_positive_numeric_datum_value":pollen["total_positive_numeric_datum_value"],
                    "this_record_accepted_as_21ka_population_occupancy":False,
                    **item})
            row["source_labelled_18_24ka_cal_bp_age_model_n"]=candidate_count
            if candidate_count:s_counts["positive_datum_samples_with_source_labelled_cal_BP_18_24ka_n"]+=1
            source_samples_rows.append(row)
        dsummary={"datasetid":did,"siteid":site[did],
                  "full_dataset_samples_total_n":len(samples),
                  "taxon1299_source_datum_age_rows_total_n":used_ages,
                  "explicit_sample_age_chronologyids":";".join(sorted(chronology_ids,key=int)),
                  "source_dataset_sha256":byid[did],
                  **dict(s_counts)}
        summaries.append(dsummary)
        totals.update(s_counts)
    assert len(summaries)==6
    # The 27 v13 occurrence candidates cannot be joined by occurrence IDs to
    # these full dataset source sampleids: these are a separate source study.
    result={"schema":"chocho_original_Neotoma_6_Alnus_pollen_datum_sample_ages_v18",
        "status":STATUS,
        "source_original_datasets":6,"source_taxonid":TAXON_ID,
        "original_v13_metadata_candidate_occurrences_unjoined":26,
        "dataset_counts":summaries,
        "source_positive_datum_samples_total":totals["positive_datum_sample_n"],
        "source_positive_datum_samples_with_labelled_cal_BP_18_24ka_n":
            totals["positive_datum_samples_with_source_labelled_cal_BP_18_24ka_n"],
        "source_labelled_candidate_age_rows_n":totals["labelled_candidate_age_model_rows_n"],
        "source_labelled_candidate_age_rows_missing_both_bounds_n":totals["candidate_missing_both_age_bounds_n"],
        "source_sample_rows_with_Alnus_taxon1299":len(source_samples_rows),
        "source_age_model_rows_for_Alnus_labeled_samples":len(source_ages_rows),
        "original_sample_ages_independently_calibrated":0,
        "source_Alnus_pollen_morphology_independently_validated":0,
        "accepted_21ka_Alnus_host_occupied_population_sites":0,
        "original_frozen_known_host_specificity_negative_results_unchanged":True,
        "source_sha256":byid,
        "limits":[
           "Direct full dataset source datum and source labelled calibrated BP models are archival metadata, not independently verified ages or species-specific pollen morphology.",
           "Multiple source age models per sample are retained; a counted source-labelled candidate sample is still one sample, not independent dates.",
           "No join to original v13 occurrence ID was verified in v17; do not impute an occurrence crosswalk by matched depths/ages.",
           "Missing age uncertainty, chronology calibration, sediment provenance and non-local transported pollen all preclude ancient occupied plant claims.",
           "The original insect genetic response and host-climate alternative real plant negative controls remain unopened/unchanged."]}
    return source_samples_rows,source_ages_rows,summaries,result

def csvwrite(path,rows):
    if not rows:raise ValueError("cannot materialize empty source data without classification")
    header=sorted(set().union(*(r.keys() for r in rows)))
    with Path(path).open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=header,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contract",required=True,type=Path)
    p.add_argument("--full-source-dir",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    c=json.loads(a.contract.read_text())
    samples,ages,datasets,summary=validate_source(c,a.full_source_dir)
    csvwrite(a.out/"all_source_Alnus_sample_datum_v18.csv",samples)
    csvwrite(a.out/"all_source_Alnus_sample_age_models_v18.csv",ages)
    csvwrite(a.out/"six_source_dataset_pollen_age_inventory_v18.csv",datasets)
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()

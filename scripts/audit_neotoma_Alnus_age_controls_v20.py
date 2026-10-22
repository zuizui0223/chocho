#!/usr/bin/env python3
"""Audit direct Neotoma chronology reliable age span and age-control depth.

Compare 21 originally source-labelled positive Alnus pollen sample ages to
their EXACT original chronologyid. Preserve all multiple-age-model variants.
Neither database model metadata nor pollen taxon label independently verifies
21ka occupied local Alnus glutinosa trees, larval hosts, or genetic breaks.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

SOURCE_SHA={
"v18_ages":"85cb546d43e57c8ff195c09d5dba253e6e2cc88a5cbae4e1eddb762e5c67c6e4",
"v18_samples":"22e2ceae388c704dcc1dc21d6e71605f674807ec3ab3dce004956b9220d674f4",
"v18_summary":"934c29c61a9ee8603f54758fcf27f85f1970bd512e39cebf9a44fdebd9a195cf",
"v19_summary":"bc70337b6cff34fac59823bd4861bb3ebed315099915bbc94e996cfc0ae67166",
"v19_routes":"d29c8c1848316d11a58014ec3a452f65b6fdfe582cd5df041d4b4e7b8ba5baa5"
}
STATUS="SOURCE_CHRONOLOGY_DEPTH_AND_RELIABILITY_CAVEATS_NOT_VERIFIED_FOSSIL_SPECIES"
CAL="calibrated radiocarbon years bp"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(1048576),b""):h.update(x)
    return h.hexdigest()

def rows(path):
    with Path(path).open("r",newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def number(value):
    if value is None or isinstance(value,bool) or str(value).strip()=="":
        return None
    try:x=float(value)
    except (ValueError,TypeError,OverflowError):return None
    return x if math.isfinite(x) else None

def chronology_object(data,expected_id,expected_dataset):
    if not isinstance(data,dict) or not isinstance(data.get("data"),list) or not data["data"]:
        raise ValueError("official chronology must include at least one source object")
    records=data["data"]
    c=records[0].get("chronology") if isinstance(records[0],dict) else None
    if not isinstance(c,dict) or str(c.get("chronologyid"))!=str(expected_id):
        raise ValueError("source chronology identity not explicitly verified")
    # Neotoma v2 sometimes returns the SAME chronology metadata in multiple
    # data[] records, with different age-control objects in each record.
    # Looking only at data[0] silently drops dated age controls and falsely
    # labels original pollen samples as beyond the known depth range.
    control_records=[]
    for entry in records:
        detail=entry.get("chronology") if isinstance(entry,dict) else None
        if not isinstance(detail,dict) or str(detail.get("chronologyid"))!=str(expected_id):
            raise ValueError("mixed chronology identities in one original official response")
        for key in ("agetype","default","modelType","reliableagespan"):
            if detail.get(key)!=c.get(key):
                raise ValueError("inconsistent repeated source chronology metadata: "+key)
        current=detail.get("controls")
        if not isinstance(current,list):
            raise ValueError("a repeated chronology lacks its own controls array")
        control_records.extend(current)
    linked=c.get("datasets")
    if not isinstance(linked,list) or not any(
            isinstance(v,dict) and str(v.get("datasetid"))==str(expected_dataset) and
            str(v.get("datasettype","")).lower()=="pollen" for v in linked):
        raise ValueError("source chronology not explicitly linked to original pollen dataset")
    age_type=str(c.get("agetype") or "")
    reliable=c.get("reliableagespan")
    if not isinstance(reliable,dict):
        raise ValueError("official chronology reliable age span absent")
    younger=number(reliable.get("younger"));older=number(reliable.get("older"))
    if (younger is not None and older is not None and younger>older):
        raise ValueError("inverted reliable age bounds")
    controls=[]
    dedup=set()
    for v in control_records:
        if not isinstance(v,dict):raise ValueError("malformed chronology control")
        key=(str(v.get("chroncontrolid")),number(v.get("depth")),number(v.get("age")))
        if key not in dedup:
            controls.append(v)
            dedup.add(key)
    cc=[]
    for v in controls:
        cc.append({"chroncontrolid":str(v.get("chroncontrolid") or ""),
                   "depth_cm":number(v.get("depth")),
                   "age_BP_numeric":number(v.get("age")),
                   "control_type":str(v.get("controltype") or ""),
                   "geochron_lab_number":((v.get("geochron") or {}).get("geochronlabnumber")
                       if isinstance(v.get("geochron"),dict) else None),
                   "dated_material":((v.get("geochron") or {}).get("materialdated")
                       if isinstance(v.get("geochron"),dict) else None),
                   "geochron_age_scale":((v.get("geochron") or {}).get("geochronagetype")
                       if isinstance(v.get("geochron"),dict) else None),
                   "geochron_raw_age_BP":number(((v.get("geochron") or {}).get("geochronage") or {}).get("age"))
                       if isinstance(v.get("geochron"),dict) else None,
                   "source_limits":v.get("controllimits")})
    validdepth=sorted([x["depth_cm"] for x in cc if x["depth_cm"] is not None])
    num_control_ages=sum(x["age_BP_numeric"] is not None for x in cc)
    return {"chronologyid":str(expected_id),"source_datasetid":str(expected_dataset),
            "agetype":age_type,
            "isdefault":c.get("default") is True,
            "model_type":str(c.get("modelType") or ""),
            "chronology_name":str(c.get("chronologyName") or ""),
            "date_prepared":c.get("datePrepared"),
            "chronology_notes":str(c.get("chronologynotes") or ""),
            "source_reliable_younger_BP":younger,"source_reliable_older_BP":older,
            "source_age_controls_n":len(controls),
            "official_chronology_response_records_n":len(records),
            "official_chronology_controls_before_dedup_n":len(control_records),
            "numeric_age_controls_n":num_control_ages,
            "controls_with_valid_depth_n":len(validdepth),
            "source_control_depth_min_cm":min(validdepth) if validdepth else None,
            "source_control_depth_max_cm":max(validdepth) if validdepth else None,
            "distinct_control_depths_n":len(set(validdepth)),
            "control_n_with_original_14C_lab_number":sum(bool(x["geochron_lab_number"]) for x in cc),
            "controls":cc}

def risk_at_sample(sample_depth,age,metadata):
    depth=number(sample_depth);year=number(age)
    if depth is None or year is None:
        raise ValueError("original sample lacks numeric source depth/central age")
    young=metadata["source_reliable_younger_BP"];old=metadata["source_reliable_older_BP"]
    within=(young<=year<=old) if young is not None and old is not None else None
    lo=metadata["source_control_depth_min_cm"]
    hi=metadata["source_control_depth_max_cm"]
    if lo is None or hi is None:depth_status="SOURCE_CONTROLS_WITHOUT_NUMERIC_DEPTH"
    elif depth<lo:depth_status="SHALLOWER_THAN_SHALLOWEST_RETURNED_CONTROL"
    elif depth>hi:depth_status="DEEPER_THAN_DEEPEST_RETURNED_CONTROL"
    else:depth_status="WITHIN_RETURNED_CONTROL_DEPTH_EXTENT"
    before=sorted((x for x in metadata["controls"] if x["depth_cm"] is not None
                   and x["age_BP_numeric"] is not None and x["depth_cm"]<=depth),
                  key=lambda r:r["depth_cm"])
    after=sorted((x for x in metadata["controls"] if x["depth_cm"] is not None
                  and x["age_BP_numeric"] is not None and x["depth_cm"]>=depth),
                  key=lambda r:r["depth_cm"])
    below=before[-1] if before else None
    above=after[0] if after else None
    return {"sample_depth_cm":depth,"sample_center_age_BP":year,
            "source_model_age_span_older_BP":old,
            "source_model_age_span_younger_BP":young,
            "sample_center_inside_database_reliable_age_span":within,
            "relative_to_source_returned_control_depth_extent":depth_status,
            "nearest_shallower_control_depth_cm":below["depth_cm"] if below else None,
            "nearest_shallower_control_age_BP":below["age_BP_numeric"] if below else None,
            "nearest_deeper_control_depth_cm":above["depth_cm"] if above else None,
            "nearest_deeper_control_age_BP":above["age_BP_numeric"] if above else None,
            "nearest_bracketing_control_age_gap_BP":abs(above["age_BP_numeric"]-below["age_BP_numeric"])
               if below and above and below["depth_cm"]<above["depth_cm"] else None,
            "exact_control_at_sample_depth":bool((below and below["depth_cm"]==depth) or
                                                 (above and above["depth_cm"]==depth))}

def flatten_controls(metadata):
    rows=[]
    for chron in metadata.values():
        for c in chron["controls"]:
            rows.append({"datasetid":chron["source_datasetid"],"chronologyid":chron["chronologyid"],
                        "source_agetypelabel":chron["agetype"],
                        "default_source_chronology":chron["isdefault"],**c})
    return rows

def analyze(contract,source_dir,ages,samples,v18,v19,receipt):
    if (contract["status"]!="EXPLORATORY_POST_RETRIEVAL_SOURCE_METADATA_QC_NOT_INDEPENDENT_CALIBRATION"
        or len(ages)!=788 or len(samples)!=437
        or v18["source_positive_datum_samples_with_labelled_cal_BP_18_24ka_n"]!=21
        or v19["chronology_endpoints_identifiable"]!=11
        or v19["fake_chronology_negative_control_is_safe"] is not True):
        raise ValueError("original v18/v19 source population ages or official chronology gate drift")
    expected=contract["sources"]["chronology_ids"]
    by_dataset={str(int(k)):str(int(k)) for k in contract["sources"]["dataset_ids"]}
    route_sha=receipt["source_bodies"]
    metadata={}
    for id_ in expected:
        filename=f"chronology_{id_}.json.gz"
        path=Path(source_dir)/filename
        if filename not in route_sha or sha(path)!=route_sha[filename]:
            raise ValueError("original v19 chronology source body SHA mismatch "+filename)
        raw=json.loads(gzip.decompress(path.read_bytes()))
        # Source dataset may also include associated geochronological project;
        # the requested fossil pollen dataset must match the v18 source.
        original_dataset=None
        for d in raw.get("data",[]):
            ch=d.get("chronology",{})
            if str(ch.get("chronologyid"))==str(id_):
                possible={str(r.get("datasetid")):r for r in ch.get("datasets",[])
                          if isinstance(r,dict) and str(r.get("datasettype","")).lower()=="pollen"}
                chosen=list(set(possible)&set(by_dataset))
                if len(chosen)!=1:raise ValueError("unexpected pollen source dataset links")
                original_dataset=chosen[0]
        if original_dataset is None:raise ValueError("chronology identity not in original linked datasets")
        metadata[str(id_)]=chronology_object(raw,id_,original_dataset)
    sample_map={(r["datasetid"],r["sampleid"]):r for r in samples}
    if len(sample_map)!=437:raise ValueError("original source pollen samples not unique")
    age_by_sample=defaultdict(list)
    for r in ages:age_by_sample[(r["datasetid"],r["sampleid"])].append(r)
    targets=[r for r in ages
        if r["source_labelled_18_24ka_cal_BP_candidate"]=="True" and r["datum_positive"]=="True"]
    if len(targets)!=21 or len({(r["datasetid"],r["sampleid"]) for r in targets})!=21:
        raise ValueError("original 21 candidate source chronology sample identities drift")
    report=[]
    for row in targets:
        key=(row["datasetid"],row["sampleid"])
        source_sample=sample_map.get(key)
        if not source_sample:raise ValueError("sample-to-age source join missing explicitly linked sampleid")
        m=metadata[row["source_chronologyid"]]
        if m["source_datasetid"]!=row["datasetid"] or m["agetype"].casefold()!=CAL:
            raise ValueError("original candidate source chronology age scale disagrees with full object")
        info=risk_at_sample(source_sample["sample_depth"],row["source_age_BP_numeric"],m)
        other=[]
        for alt in age_by_sample[key]:
            if alt is row:continue
            chid=alt["source_chronologyid"]
            other.append({"chronologyid":chid,
                "original_age_numeric":number(alt.get("source_age_BP_numeric")),
                "original_age_type":alt.get("source_agetype_label"),
                "default_model_flag_from_chronology":metadata[chid]["isdefault"] if chid in metadata else None})
        if not any(alt["source_chronologyid"]==row["source_chronologyid"] for alt in age_by_sample[key]):
            raise ValueError("candidate original age lost")
        report.append({"datasetid":row["datasetid"],"siteid":row["siteid"],
            "sampleid":row["sampleid"],"original_Alnus_pollen_positive_datum_value":
                number(row["source_positive_numeric_datum_value"]),
            "source_chronologyid":row["source_chronologyid"],
            "source_age_type":row["source_agetype_label"],
            "source_age_older_bound_BP":number(row.get("source_ageolder_BP")),
            "source_age_younger_bound_BP":number(row.get("source_ageyounger_BP")),
            "source_sample_has_both_age_bounds":row["has_both_source_age_bounds"]=="True",
            "original_chronology_isdefault":m["isdefault"],
            "original_chronology_model_type":m["model_type"],
            "chronology_control_n":m["source_age_controls_n"],
            "chronology_dated_lab_controls_n":m["control_n_with_original_14C_lab_number"],
            "alternative_original_source_age_models":other,
            "only_source_metadata_not_independent_species_occupancy":True,
            **info})
    if len(report)!=21:raise ValueError("not 21 exact source-sample records")
    counts=Counter()
    for r in report:
        counts["source_labelled_samples"]+=1
        counts["nondefault_source_chronology"]+=not r["original_chronology_isdefault"]
        counts["source_sample_missing_both_age_bounds"]+=not r["source_sample_has_both_age_bounds"]
        counts["beyond_oldest_reliable_age_span"]+=r["sample_center_inside_database_reliable_age_span"] is False
        counts["in_database_reliable_age_span"]+=r["sample_center_inside_database_reliable_age_span"] is True
        counts["deeper_than_last_returned_source_control"]+=r["relative_to_source_returned_control_depth_extent"]=="DEEPER_THAN_DEEPEST_RETURNED_CONTROL"
        counts["within_returned_control_depth_extent"]+=r["relative_to_source_returned_control_depth_extent"]=="WITHIN_RETURNED_CONTROL_DEPTH_EXTENT"
        counts["exact_dated_control_at_sample_depth"]+=r["exact_control_at_sample_depth"]
        counts["bracket_control_age_gap_over_10000_BP"]+=(
            r["nearest_bracketing_control_age_gap_BP"] is not None and
            r["nearest_bracketing_control_age_gap_BP"]>10000)
    by_site=defaultdict(list)
    for r in report:by_site[r["datasetid"]].append(r)
    sites=[{"datasetid":d,"n":len(x),
       "nondefault_model_samples":sum(not r["original_chronology_isdefault"] for r in x),
       "outside_reported_reliable_age_span":sum(r["sample_center_inside_database_reliable_age_span"] is False for r in x),
       "deeper_than_deepest_source_control":sum(r["relative_to_source_returned_control_depth_extent"]=="DEEPER_THAN_DEEPEST_RETURNED_CONTROL" for r in x),
       "bracketed_in_returned_control_depth_extent":sum(r["relative_to_source_returned_control_depth_extent"]=="WITHIN_RETURNED_CONTROL_DEPTH_EXTENT" for r in x)}
       for d,x in sorted(by_site.items(),key=lambda t:int(t[0]))]
    overview={"schema":"chocho_Neotoma_Alnus_21ka_sample_depth_vs_direct_chroncontrols_v20",
       "status":STATUS,"original_v18_source_labelled_18_24ka_positive_sample_n":21,
       "original_v19_official_chronologies":len(metadata),
       "chronology_metadata_is_not_independent_calibration":True,
       "original_sample_ages_independently_calibrated":0,
       "source_Alnus_glutinosa_pollen_species_independently_validated":0,
       "ancient_21ka_Alnus_occupied_plant_populations_accepted":0,
       "genetic_response_used":False,
       "descriptive_counts":dict(counts),"dataset_status":sites,
       "chronology_ids":sorted(metadata,key=int),
       "source_sha256":dict(SOURCE_SHA),
       "limits":["18ka-24ka source labels are age model predictions, not independent 14C calibration or direct dated specimens",
         "Depth beyond deepest returned control is evidence of possible extrapolation, not proof the external age model is wrong",
         "Inside source reliableage span is not a substitute for missing individual 95% age bounds",
         "Multiple model chronologies are not combined, arbitrarily shifted, or selected for positive LGM evidence",
         "Direct Alnus taxonid1299 pollen data are not an independently species-diagnostic fossil or occupied plant population",
         "No historic plant resource connectivity or insect genetic barrier established, and earlier negative real-other-plant TTF/chocho controls remain unchanged"]}
    return list(metadata.values()),flatten_controls(metadata),report,sites,overview

def csvwrite(path,records):
    if not records:raise ValueError("empty source record list")
    fields=sorted(set().union(*(r.keys() for r in records)))
    with Path(path).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        for r in records:
            w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v
                        for k,v in r.items()})

def main():
    p=argparse.ArgumentParser()
    for name in ("contract","source-dir","v18-ages","v18-samples",
                  "v18-summary","v19-summary","v19-receipt","v19-routes","out"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    for key,path in (("v18_ages",a.v18_ages),("v18_samples",a.v18_samples),
                     ("v18_summary",a.v18_summary),("v19_summary",a.v19_summary),
                     ("v19_routes",a.v19_routes)):
        if sha(path)!=SOURCE_SHA[key]:raise ValueError("frozen v18/v19 source SHA changed "+key)
    contract=json.loads(a.contract.read_text())
    receipt=json.loads(a.v19_receipt.read_text())
    if len(receipt["source_bodies"])<17:raise ValueError("v19 official chronology payloads missing")
    a.out.mkdir(parents=True,exist_ok=True)
    models,controls,candidates,sites,summary=analyze(
        contract,a.source_dir,rows(a.v18_ages),rows(a.v18_samples),
        json.loads(a.v18_summary.read_text()),json.loads(a.v19_summary.read_text()),receipt)
    csvwrite(a.out/"all_11_original_chronology_metadata_v20.csv",models)
    csvwrite(a.out/"all_source_chronology_controls_v20.csv",controls)
    csvwrite(a.out/"all_21_original_Alnus_age_depth_caveats_v20.csv",candidates)
    csvwrite(a.out/"five_source_dataset_age_support_v20.csv",sites)
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()

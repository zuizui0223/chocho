#!/usr/bin/env python3
"""Audit six source-bound Neotoma Alnus 18–24ka pollen datasets for age SCALE.

A source occurrence's numeric age can disagree with its dataset's published
age range, and Neotoma holds mixed radiocarbon/calibrated-year chronologies.
That is a source-identity/chronology review FLAG, not proof a fossil is false.

Only original already retrieved six-dataset metadata and 27 candidate source
records are accessed. No new plant/fossil/genetic data, no 21ka occupancy.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path

SOURCES={
  "candidate":"3d47fe9a31cb6a04ce60aa0f1d7fc8dd7fb1508d3966d6f3a7dc6942ddfa9efc",
  "candidate_summary":"16bc59d362365416ac9144d29bae05ed503341ebb99d109f750ababdaa401bd9",
  "datasets":"c212dd7a3e8cb888097a47575bb3c7c2007add5735c96d47c7f5c33888e8ecef",
  "dataset_summary":"8000b0a79a9afb93740b0698db4b2703dd5b93d749373777bb2dd50add77696b"
}
EXPECTED_DATASET_SITES={
  "4396":"3392","14556":"9834","22662":"14474",
  "22683":"14485","24197":"15750","41351":"23855"
}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda:f.read(1024*1024),b""):h.update(buf)
    return h.hexdigest()

def csvrows(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def source_previews(dataset_rows):
    """Read only source previews; chronology object may be truncated at 5,500 chars."""
    groups=defaultdict(dict)
    for r in dataset_rows:
        did=r["datasetid"]
        if did=="9999999":continue
        if (did not in EXPECTED_DATASET_SITES or
            r["siteid"]!=EXPECTED_DATASET_SITES[did] or
            r["status"]!="SOURCE_METADATA_AVAILABLE_NOT_INDEPENDENTLY_VALIDATED"):
            raise ValueError("original six dataset/source IDs or access statuses changed")
        typ=r["source_type"]
        if typ not in ("dataset","chronologies","publications") or typ in groups[did]:
            raise ValueError("missing, duplicate or unexpected source endpoint")
        groups[did][typ]=r
    if len(groups)!=6 or any(len(x)!=3 for x in groups.values()):
        raise ValueError("not the exact six official dataset/chronology/publication source triples")
    return groups

def _preview_string(source,text_key):
    # FIRST chronology only, source preview may be truncated; no later chronology assertions.
    # ChronologyName etc are plain source quoted strings, possibly \\" escapes.
    matched=re.search(r'"'+re.escape(text_key)+r'":\s*("(?:[^"\\]|\\.)*"|null)',source)
    if not matched:return None
    return None if matched.group(1)=="null" else json.loads(matched.group(1))

def parse_dataset_source(did,source_triple):
    item=source_triple["dataset"]
    raw=json.loads(item["first_record_preview"])
    site=raw.get("site")
    if not isinstance(site,dict):raise ValueError("source dataset preview missing original site")
    dataset=[a for a in site.get("datasets",[]) if str(a.get("datasetid"))==did]
    if len(dataset)!=1 or str(site.get("siteid"))!=EXPECTED_DATASET_SITES[did]:
        raise ValueError("Neotoma dataset and site source metadata crosswalk drift")
    d=dataset[0]
    if d.get("datasettype")!="pollen":raise ValueError("original source proxy no longer pollen")
    ranges=d.get("agerange")
    if not isinstance(ranges,list) or not ranges:
        raise ValueError("no documented source dataset age ranges")
    range_values=[]
    for r in ranges:
        start,end=r.get("ageyoung"),r.get("ageold")
        unit=r.get("units")
        if (not isinstance(start,(int,float)) or not isinstance(end,(int,float))
            or end<start or not isinstance(unit,str) or not unit.strip()):
            raise ValueError("uninterpretable original source chronology unit or numeric range")
        range_values.append({"start":start,"end":end,"unit":unit})
    source_chrono=source_triple["chronologies"]
    chronology_preview=source_chrono["first_record_preview"]
    chrono_agetype=_preview_string(chronology_preview,"agetype")
    chrono_name=_preview_string(chronology_preview,"chronologyName")
    chrono_note=_preview_string(chronology_preview,"chronologynotes")
    publication=source_triple["publications"]["first_record_preview"]
    p=json.loads(publication).get("publication")
    if not isinstance(p,dict) or not p.get("citation"):
        raise ValueError("source publication citation unavailable")
    # Age range units and first chronology are heterogeneous; do not choose a chronology.
    units={r["unit"] for r in range_values}
    flags=[]
    if len(units)>1:
        flags.append("MULTIPLE_INCOMPATIBLE_PUBLISHED_AGE_SCALES")
    if chrono_agetype and chrono_agetype not in units:
        flags.append("FIRST_CHRONOLOGY_AGE_TYPE_NOT_IN_DATASET_RANGE_UNITS")
    if chrono_note and "provisional" in chrono_note.lower():
        flags.append("SOURCE_SAYS_CHRONOLOGY_VERY_PROVISIONAL")
    if chrono_note and ("IntCal04" in chrono_note or "intcal04" in chrono_note.lower()):
        flags.append("OLD_CALIBRATION_CURVE_SOURCE_REVIEW")
    return {"datasetid":did,"siteid":EXPECTED_DATASET_SITES[did],
      "site_name":str(site.get("sitename") or ""),
      "dataset_type":d["datasettype"],"dataset_age_ranges":range_values,
      "dataset_age_unit_labels":";".join(sorted(units)),
      "dataset_contains_calibrated_BP":any("calibrated" in x.lower() for x in units),
      "dataset_contains_uncalibrated_radiocarbon_BP":any(
           "radiocarbon years BP" in x.lower() and "calibrated" not in x.lower()
           for x in units),
      "first_chronology_age_type":chrono_agetype,
      "first_chronology_name":chrono_name,
      "first_chronology_notes":chrono_note,
      "number_of_source_chronology_entries":int(source_chrono["data_count"]),
      "publication_citation":p["citation"],
      "dataset_doi":";".join(d.get("doi") or []),
      "metadata_flags":flags}

def align_records(candidate_rows,by_dataset):
    if len(candidate_rows)!=27:
        raise ValueError("not original 27 Alnus candidate source rows")
    res=[];ids=set();per_dataset=defaultdict(list)
    for r in candidate_rows:
        did=r["source_datasetid"]
        if (did not in by_dataset or r["source_siteid"]!=by_dataset[did]["siteid"]
            or r["source_sample_taxonid"]!="1299" or
            r["source_sample_taxonname"]!="Alnus glutinosa" or
            r["source_datasettype"]!="pollen" or
            r["source_grade"]!="SOURCE_18_24KA_METADATA_CANDIDATE_NOT_VERIFIED_FOSSIL_OCCUPANCY"):
            raise ValueError("original fossil source proxy, site, taxon ID or candidate grade drift")
        age=float(r["sample_age_BP"])
        if not 18000<=age<=24000:
            raise ValueError("prior exact numeric source age was not 18–24ka")
        ranges=by_dataset[did]["dataset_age_ranges"]
        in_any=any(x["start"]<=age<=x["end"] for x in ranges)
        older=r["sample_age_older_BP"]
        younger=r["sample_age_younger_BP"]
        status=("AGE_VALUE_OUTSIDE_ALL_PUBLISHED_DATASET_AGE_ENVELOPES_REVIEW_UNITS"
                if not in_any else "WITHIN_AT_LEAST_ONE_DATASET_ENVELOPE_NOT_CHRONOLOGY_VERIFIED")
        flags=list(by_dataset[did]["metadata_flags"])
        if not in_any:flags.append("SOURCE_AGE_OUTSIDE_DATASET_ADVERTISED_ENVELOPE")
        if not (older and younger):flags.append("SAMPLE_AGE_UNCERTAINTY_NOT_IN_OCCURRENCE_RESPONSE")
        ident=r["occid"]
        if not ident:raise ValueError("source fossil candidate missing occurrence ID")
        if ident in ids:flags.append("DUPLICATE_OCCURRENCE_NOT_INDEPENDENT_OBSERVATION")
        ids.add(ident)
        entry={"datasetid":did,"siteid":r["source_siteid"],
          "site_name":by_dataset[did]["site_name"],
          "occid":ident,"sample_age_BP_as_returned":age,
          "source_age_status":status,
          "source_age_range_envelope_compatibility_numeric_only":in_any,
          "sample_age_younger_BP_as_returned":younger or None,
          "sample_age_older_BP_as_returned":older or None,
          "first_chronology_age_type":by_dataset[did]["first_chronology_age_type"],
          "dataset_age_unit_labels":by_dataset[did]["dataset_age_unit_labels"],
          "flags":flags,
          "fossil_species_occurrence_validated":False,
          "sample_age_calibration_model_verified":False}
        res.append(entry);per_dataset[did].append(entry)
    if len(ids)!=26:
        raise ValueError("original 26 distinct Neotoma occurrence IDs changed")
    summary={}
    for did,rows in per_dataset.items():
        out=sum(x["source_age_status"].startswith("AGE_VALUE_OUTSIDE") for x in rows)
        summary[did]={"candidate_rows":len(rows),
                      "distinct_occurrence_ids":len({x["occid"] for x in rows}),
                      "numeric_age_outside_any_published_range":out,
                      "numeric_age_range_min":min(x["sample_age_BP_as_returned"] for x in rows),
                      "numeric_age_range_max":max(x["sample_age_BP_as_returned"] for x in rows)}
    if set(summary)!=set(EXPECTED_DATASET_SITES):
        raise ValueError("some of six original source datasets no longer represented")
    return res,summary

def run(args):
    for k,source in (("candidate",args.candidate),("candidate_summary",args.candidate_summary),
                     ("datasets",args.dataset_metadata),("dataset_summary",args.dataset_summary)):
        if sha(source)!=SOURCES[k]:raise ValueError("original source SHA changed: "+k)
    s=json.loads(args.candidate_summary.read_text())
    meta=json.loads(args.dataset_summary.read_text())
    if (s["source_candidate_distinct_occids"]!=26 or
        meta["source_query_count"]!=19 or
        meta["candidate_source_datasets"]!=6 or
        meta["sample_ages_independently_chronology_validated"] or
        meta["pollen_morphology_species_level_validated"]):
        raise ValueError("original source chronology uncertainty or sample cohort drift")
    sources=source_previews(csvrows(args.dataset_metadata))
    by={did:parse_dataset_source(did,v) for did,v in sources.items()}
    rows,counts=align_records(csvrows(args.candidate),by)
    out_n=sum(x["numeric_age_outside_any_published_range"] for x in counts.values())
    flags=Counter(f for x in by.values() for f in x["metadata_flags"])
    source_grade={
      "schema":"chocho_six_Alnus_18_24ka_candidate_pollen_source_age_units_QC_v15",
      "status":"SOURCE_AGE_SCALE_CONFLICTS_REQUIRING_FOSSIL_CHRONOLOGY_AND_TAXONOMY_REVIEW",
      "original_dataset_count":6,"original_source_candidate_rows":27,
      "unique_candidate_occurrence_ids":26,
      "candidate_rows_outside_all_published_dataset_age_envelopes":out_n,
      "by_dataset_id":counts,
      "site_metadata_risk_flags":dict(flags),
      "every_dataset_chronology_or_publication_accessible":True,
      "source_sample_calibrated_18_24ka_ages_verified_n":0,
      "independent_species_pollen_identifications_verified_n":0,
      "fossil_occurrences_admitted_to_LGM_host_occupancy_model_n":0,
      "insect_genetic_response_opened":False,
      "original_true_host_climate_specificity_negative_controls_preserved":True,
      "source_sha256":dict(SOURCES),
      "limits":["Age values outside coarse dataset age ranges are flags, not definitive invalid fossils: chronologies and calibrated/un-calibrated units may differ",
        "The metadata preview covers ONLY first source chronology, sometimes truncated; does not identify which chronology is assigned to the sample",
        "Presence of a chronology/control age does not imply candidate pollen sample age posterior is independently verified",
        "Mixed radiocarbon versus calibrated years BP cannot be compared directly or converted by subtraction of a fixed offset",
        "Even exact Neotoma binomial pollen labels do not establish species diagnostic morphology or ancient occupied populations",
        "Six technical-control datasets are not independent fossil evidence for 416 plant species or LGM insect genetic barriers"]}
    if out_n!=4 or counts["14556"]["numeric_age_outside_any_published_range"]!=3 or counts["22662"]["numeric_age_outside_any_published_range"]!=1:
        raise ValueError("source original four coarse age-envelope mismatches were not reproduced")
    return rows,list(by.values()),source_grade

def write_csv(p,rows):
    if not rows:raise ValueError("empty audit")
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with Path(p).open("w",newline="",encoding="utf-8") as stream:
        w=csv.DictWriter(stream,fieldnames=keys,lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(v,sort_keys=True,ensure_ascii=False) if isinstance(v,(dict,list))
                        else v for k,v in r.items()})

def main():
    p=argparse.ArgumentParser()
    for k in ("candidate","candidate-summary","dataset-metadata","dataset-summary","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    r,s,q=run(a)
    write_csv(a.out/"all_27_candidate_source_age_scale_flags_v15.csv",r)
    write_csv(a.out/"six_dataset_source_chronology_risk_v15.csv",s)
    (a.out/"summary.json").write_text(json.dumps(q,indent=2,sort_keys=True)+"\n")
    print(json.dumps(q,sort_keys=True))
if __name__=="__main__":main()

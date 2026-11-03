#!/usr/bin/env python3
"""Version the exact original response-blind 500 GBIF indexed plant×island
queries after source-preserving fixes of rate-limit errors and rejected
multi-island WKT. No rows added, no missing results silently imputed, and
no historically indexed year treated as the biological plant introduction.
"""
from __future__ import annotations
import argparse,csv,json,statistics
from collections import Counter
from pathlib import Path

def read(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def save(p,rows):
    if not rows:raise RuntimeError("No source rows")
    fields=[]
    for r in rows:
        for field in r:
            if field not in fields:fields.append(field)
    with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(rows)

def key(row,code="island_code"):
    return (str(row["plant_id"]),str(row[code]))

def main():
    p=argparse.ArgumentParser()
    for x in ("original500","original-audit","429-retry","initial-split","remaining-split","out-dir"):
        p.add_argument("--"+x,type=Path,required=True)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    rows=read(a.original500)
    original=json.loads(a.original_audit.read_text(encoding="utf-8"))
    retry=read(a.__dict__["429_retry"])
    first=read(a.initial_split)
    rest=read(a.remaining_split)
    if len(rows)!=500 or original["pilot_dated_pairs"]!=219 or len(retry)!=12:
        raise RuntimeError("Exact source-fixed 500-pair GBIF cohort changed")
    if len({key(r) for r in rows})!=500:
        raise RuntimeError("Duplicate original host×WGSRPD3 key")
    originals={key(r):r for r in rows}
    if sum(r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND" for r in rows)!=219:
        raise RuntimeError("Original 219 dated source status drift")
    extras={}
    def update(row,status,year,records,source):
        k=key(row)
        if k not in originals or originals[k]["q4_pilot_status"]!="GBIF_QUERY_FAILED":
            raise RuntimeError(f"Rescue must only update original failed status {k}")
        if status not in {"DATED_OCCURRENCE_FOUND","NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"}:
            raise RuntimeError("Not a successful source outcome")
        if status=="DATED_OCCURRENCE_FOUND" and not year:
            raise RuntimeError("Dated plant rescue lacks earliest source indexed year")
        if status=="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT" and (year or str(records)!="0"):
            raise RuntimeError("No-index status with nonzero records")
        if k in extras:raise RuntimeError("Multiple rescues for same fixed source cell "+str(k))
        extras[k]={"q4_pilot_status":status,
                   "earliest_gbif_indexed_year":str(year),
                   "gbif_plant_record_count":str(records),
                   "source_observation_repair_stage":source}
    for r in retry:
        if r["status"] in {"DATED_OCCURRENCE_FOUND","NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"}:
            update(r,r["status"],r.get("earliest_indexed_year") or "",r.get("n_records") or "",
                   "ONE_POLYGON_OR_SMALL_MULTI_ORIGINAL_HTTP429_RETRY")
    for stage, arr in [
        ("VALIDATED_FULL_POLYGON_INITIAL_HTTP400",first),
        ("VALIDATED_FULL_POLYGON_REMAINING_HTTP400",rest)
    ]:
        for r in arr:
            if r.get("kind") not in ("RESCUE","RESCUE_REMAINDER"):continue
            if str(r.get("complete_full_WGSRPD_region")).lower()!="true":continue
            update({"plant_id":r["plant_id"],"island_code":r["region"]},r["status"],
                   r.get("earliest_indexed_year_if_complete") or "",
                   r.get("total_records_if_complete") or "",stage)
    stage_counts=Counter(item["source_observation_repair_stage"] for item in extras.values())
    merged=[]
    for r in rows:
        r={**r}
        k=key(r)
        r["source_observation_repair_stage"]="ORIGINAL_FIXED_500"
        if k in extras:
            patch=extras[k]
            r.update(patch)
            # Previously failed queries are now valid source first-year values.
            # Detailed per-year counts are not available from split rescue.
            if "FULL_POLYGON" in patch["source_observation_repair_stage"]:
                for field in ("gbif_plant_pre1900_count","gbif_plant_pre1950_count",
                              "gbif_plant_pre2001_count","distinct_indexed_years",
                              "year_facet_complete"):
                    r[field]=""
            r["error"]=""
        merged.append(r)
    count=Counter(r["q4_pilot_status"] for r in merged)
    dated=[int(r["earliest_gbif_indexed_year"]) for r in merged
           if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND"]
    if len(dated)!=219+sum(p["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND" for p in extras.values()):
        raise RuntimeError("Merged source dates missing after rescue")
    if sum(count.values())!=500:raise RuntimeError("Incorrect source population")
    byisland=Counter(r["island_code"] for r in merged if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND")
    payload={
      "schema":"chocho_Q4_same_500_host_x_island_historical_GBIF_date_revision_v0.2",
      "status":"TECHNICALLY_RESCUED_INDEX_DATES_NOT_IMPORT_OR_RESIDENCE_AGE",
      "original_fixed_blind_source_pairs":500,
      "original_source_dated":219,
      "new_rescued_pairs":len(extras),
      "new_rescued_dated":sum(r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND" for r in extras.values()),
      "new_rescued_nonindexed":sum(r["q4_pilot_status"]=="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT" for r in extras.values()),
      "source_recovery_counts":dict(stage_counts),
      "original_source_status_counts":dict(Counter(r["q4_pilot_status"] for r in rows)),
      "v02_source_status_counts":dict(count),
      "v02_dated_plant_island_pairs":count["DATED_OCCURRENCE_FOUND"],
      "v02_dated_fraction_of_all_fixed_500":count["DATED_OCCURRENCE_FOUND"]/500,
      "v02_island_regions_with_dated_plant":len(byisland),
      "v02_first_indexed_year_min":min(dated),
      "v02_first_indexed_year_median":statistics.median(dated),
      "v02_first_indexed_year_max":max(dated),
      "no_new_pairs_or_butterfly_conditioned_sampling":True,
      "limits":[
         "GBIF first indexed botanical documentation is not the true first arrival or naturalization year.",
         "The first indexed year places an upper bound on introduction year and lower bound on elapsed residence time, subject to taxonomic and locality errors.",
         "Unsuccessful GBIF year searches remain technical missingness, NEVER confirmed absence.",
         "Split-polygon source requests were validated against a positive control; components cover the original official WGSRPD3 region, but exact component-boundary occurrences could double-count.",
         "The repaired sample is the same 500 response-blind plant×island opportunities; recovery fraction cannot be extrapolated as a random sample of all 5,611 or 5,545 population pairs.",
         "GBIF historical botanical effort, geography and species taxonomy produce biased documentation years; no inference of actual butterfly establishment, port transport or climate suitability."
      ]
    }
    save(a.out_dir/"fixed500_gbif_first_years_v02.csv",merged)
    (a.out_dir/"fixed500_source_recovery_audit_v02.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(payload,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()

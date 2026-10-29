#!/usr/bin/env python3
"""Pre-outcome source visits only: 1987–95 training after BSBI <=1986 host map."""
from __future__ import annotations
import argparse,collections,json
from pathlib import Path
from audit_uk_ltmn_historical_event_join_20261028 import FIELDS_VISITS,get_source,excel_rows,parse_event
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk

WINDOWS={"training_after_host_map":(1987,1995),"test":(2005,2021)}
SCHEMA="chocho_uk_LTMN_original_1987_95_host_cutoff_site_eligibility_20261029_v01"

def select(visits,cross):
    seen=set()
    per=collections.defaultdict(lambda:{x:set() for x in WINDOWS})
    years=collections.defaultdict(lambda:{x:set() for x in WINDOWS})
    for row in visits:
        key=parse_event(row)
        if key is None or key in seen:continue
        seen.add(key)
        site,year,month,day,clock=key
        if site not in cross or not 4<=month<=9:continue
        for group,(lo,hi) in WINDOWS.items():
            if lo<=year<=hi:
                per[site][group].add(key);years[site][group].add(year)
    eligible={site for site,v in per.items()
              if all(len(v[g])>=10 and len(years[site][g])>=3 for g in WINDOWS)}
    return eligible,per,years,len(seen)

def run(output):
    cross,qc=read_crosswalk()
    sites,per,years,allvisits=select(excel_rows(get_source("visits"),FIELDS_VISITS),cross)
    out={"schema":SCHEMA,
        "original_source_ukbms_visits_SHA256":"c23945feb1568b4b60ac139136c485c234699710eae395a104e0be92a958a84d",
        "original_source_ltmn_crosswalk_SHA256":"e38853b775ecadf39332a8412e2a4a2fa3c83512cf7218e4f1eda172dd83c9cd",
        "observational_windows":WINDOWS,
        "all_source_visit_events":allvisits,
        "n_geocoded_survey_transects_with_min10_3yrs_both_eras":len(sites),
        "n_visit_events_training":sum(len(per[x]["training_after_host_map"]) for x in sites),
        "n_visit_events_test":sum(len(per[x]["test"]) for x in sites),
        "n_monitoring_years_training_by_site":sorted([len(years[s]["training_after_host_map"]) for s in sites]),
        "n_monitoring_years_test_by_site":sorted([len(years[s]["test"]) for s in sites]),
        "no_original_site_ids_geocoordinates_or_butterfly_counts_exported":True,
        "outcomes_not_used_in_site_selection":True,
        "source_info_gate":"PASS_MIN10_HISTORIC_SITES" if len(sites)>=10 else "FAIL_FEWER_THAN_10_HISTORIC_SITES",
        "historic_host_data_cutoff_1986_before_all_training_visits":True,
        "status":"ORIGINAL_SOURCE_VISIT_ONLY_HISTORICAL_HOST_PREDICTION_ELIGIBILITY"}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:out[k] for k in ("status","n_geocoded_survey_transects_with_min10_3yrs_both_eras",
        "n_visit_events_training","n_visit_events_test","source_info_gate")}))
    return out

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--output",required=True,type=Path)
    run(a.parse_args().output)

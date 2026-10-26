#!/usr/bin/env python3
"""Source-duration QC after independent Henry et al. published >10-hour exclusion.

Pure data-quality readout: no niche, competition or count effect fitted.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import tempfile
from analyze_sbus_illinois_niche_counts_20261025 import frozen_focal,load_source,read_events,PROGRAMS

def run():
    focal=frozen_focal()
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/"US_butterflies_original.zip"
        load_source(path)
        events,_=read_events(path,focal)
    data={}
    for pr in PROGRAMS:
        ev=[r for (program,event),r in events.items() if program==pr]
        bad=[r for r in ev if r["bad_duration"]]
        high=[r for r in ev if not r["bad_duration"] and r["duration"]>600]
        data[pr]={
            "original_events":len(ev),
            "invalid_or_inconsistent_duration_events":len(bad),
            "duration_over_600minutes_events":len(high),
            "fraction_over_600min_of_original_events":len(high)/len(ev),
            "site_years_containing_overlong_duration":len({(r["site"],r["year"]) for r in high}),
            "duration_over_600min_threshold_is_from_Henry_2025_usage_notes":True,
            "site_names_ids_dates_and_individual_butterfly_rows_not_exported":True,
        }
    return {"schema":"chocho_sbus_effort_source_max10h_qc_20261026_v01",
            "status":"SOURCE_ONLY_DURATION_QC_NO_ECOLOGICAL_EFFECT",
            "source_dataset":"10.6084/m9.figshare.27934602.v1",
            "source_publication":"Henry et al. 2025 10.1038/s41597-025-05513-8",
            "exclusion_guideline":"Review original Edwards national trend exclusion <=0 or >10h continuous observation; no changing ecological response based on sign.",
            "programs":data}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(data,indent=2,ensure_ascii=False))

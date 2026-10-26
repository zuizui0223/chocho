#!/usr/bin/env python3
"""Check fixed Illinois host-sharing temporal signal against >10h effort outliers.

Only the 600-minute source QC changes. Original butterfly and host records,
site/year eligibility, 640 pair identities and 199 null permutations are fixed.
No site identities or survey values exported.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import tempfile

from analyze_sbus_illinois_niche_counts_20261025 import (
    frozen_focal,load_source,read_events
)
from analyze_sbus_illinois_host_degree_null_20261025 import (
    accepted_edges,butterfly_families
)
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import (
    field_panel,pair_correlations,meancontrast,null_test
)

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_RECEIPT=ROOT/"provenance/reviewer_defenses/results/sbus_illinois_adult_codetection_hostnull_20261026.json"
ORIGINAL_OBSERVED=0.14323064389607087
THRESHOLD_MIN=600


def run(hostcsv):
    first=json.loads(ORIGINAL_RECEIPT.read_text(encoding="utf-8"))
    baseline_result=first["within_site_temporal"]["observed"]
    if(abs(baseline_result["difference"]-ORIGINAL_OBSERVED)>1e-10 or
       first["sample"]["within_site_temporal_pair_count"]!=640):
        raise ValueError("Source original temporal receipt drift")
    focal=frozen_focal()
    fam=butterfly_families(focal)
    hostmap,edges,groups,viable,plantfam=accepted_edges(hostcsv,focal,fam)
    with tempfile.TemporaryDirectory() as tmp:
        raw=Path(tmp)/"source.zip"
        load_source(raw)
        events,_=read_events(raw,focal)
    baseline,baseline_keys,baseline_species,baseline_vectors,baseline_field=field_panel(events)
    original_pairs=pair_correlations(baseline,baseline_keys,baseline_species,baseline_vectors)
    if(original_pairs["status"]!="FIXED_SPECIES_PAIRS_READY"
       or len(original_pairs["within_site_temporal_candidate_pairs"])!=640):
        raise ValueError("Original field temporal roster drifted")

    overlong=[k for k,e in events.items() if k[0]=="Illinois Butterfly Monitoring Network"
              and not e["bad_duration"] and e["duration"]>THRESHOLD_MIN]
    result={
        "schema":"chocho_illinois_source_10h_effort_sensitivity_v01_20261026",
        "status":"UNVERIFIED",
        "source_data_doi":"10.6084/m9.figshare.27934602.v1",
        "source_data_original_md5":"3161b1c8802a0d67771b2044e18da928",
        "criterion":"Exclude noninvalid Illinois survey events with duration >600 min, as in published Henry 2025 source Usage Notes; original invalid duration events already excluded",
        "removed_overlong_events":len(overlong),
        "source_original_sites":baseline_field["qualified_sites"],
        "source_original_site_years":baseline_field["qualified_site_years"],
        "original_fixed_temporal_pairs":640,
        "original_temporal_shared_minus_nonshared":ORIGINAL_OBSERVED,
        "source_original_first_seed_null_q025_median_q975":
            first["within_site_temporal"]["null_difference_q025_median_q975"],
        "source_original_first_seed_right_tail":
            first["within_site_temporal"]["right_tail_fraction_exploratory"],
        "source_site_event_counts_not_persisted":True,
    }
    if len(overlong)!=5:
        result["status"]="SOURCE_DID_NOT_MATCH_PREAUDITED_DURATION_COUNT"
        return result
    for key in overlong:events.pop(key)
    try:
        revised,revised_keys,revised_species,revised_vectors,new_field=field_panel(events)
    except ValueError:
        result["status"]="NOT_IDENTICAL_COHORT_AFTER_SOURCE_QC"
        result["source_frame_fails_original_121x851_assertion"]=True
        return result
    result["qualified_after_duration_qc"]=new_field
    old_pairs=original_pairs["within_site_temporal_candidate_pairs"]
    new_pairs=pair_correlations(revised,revised_keys,revised_species,revised_vectors)
    corrected=new_pairs.get("within_site_temporal_candidate_pairs",[])
    old_set={(a,b) for a,b,r in old_pairs}
    new_set={(a,b) for a,b,r in corrected}
    result["qualified_temporal_pair_count_after_qc"]=len(corrected)
    result["original_pair_roster_unchanged"]=old_set==new_set
    result["source_eligible_site_year_effort_changed_n"]=sum(
        abs(revised[k][0]-baseline[k][0])>1e-9 for k in baseline
    )
    if(new_pairs["status"]!="FIXED_SPECIES_PAIRS_READY" or
       revised_species!=baseline_species or
       old_set!=new_set):
        result["status"]="NOT_IDENTICAL_TEMPORAL_PAIR_ROSTER_AFTER_SOURCE_QC"
        return result
    diag=meancontrast(corrected,hostmap,index=2)
    if diag["status"]!="OK" or diag["shared"]!=62 or diag["nonshared"]!=578:
        result["status"]="SOURCE_TAXON_PAIR_HOST_CLASSIFICATION_CHANGED"
        return result
    null=null_test(corrected,diag,hostmap,edges,groups,viable,plantfam,fam,2)
    result["corrected_observed"]=diag
    result["corrected_source_network_null"]=null
    result["changed_vs_original"]=diag["difference"]-ORIGINAL_OBSERVED
    result["status"]=(
        "SOURCE_10H_CORRECTION_FULL_REPLAY_QUALIFIED"
        if null["status"]=="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
        else "SOURCE_10H_CORRECTION_NULL_NOT_QUALIFIED"
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--accepted-host-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=run(a.accepted_host_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":r["status"],
        "original":r["original_temporal_shared_minus_nonshared"],
        "removed_events":r["removed_overlong_events"],
        "changed_site_year_effort":r.get("source_eligible_site_year_effort_changed_n"),
        "cohort_matches":r.get("original_pair_roster_unchanged"),
        "revised":r.get("corrected_observed"),
        "revised_null":{k:v for k,v in r.get("corrected_source_network_null",{}).items()
                        if k not in {"null_difference_history"}},
    },ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()

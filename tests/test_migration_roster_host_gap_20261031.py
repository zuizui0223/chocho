"""Verify source-locked indirect migrant-roster diagnostics without network I/O."""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031"
SPEC=importlib.util.spec_from_file_location("maritime_migrant_2025_analysis",D/"analyze_migration_roster_host_gap.py")
assert SPEC is not None and SPEC.loader is not None
MOD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def close(a,b):
    if a is None or b is None:
        assert a is None and b is None
    else:
        assert math.isclose(a,b,rel_tol=0,abs_tol=1e-9), (a,b)


def test_source_coverage_and_invalid_absence_interpretation():
    rows,source=MOD.join()
    assert len(rows)==193 and len(source)==225
    assert sum(x["source_migration"]==MOD.GROUPS[0] for x in rows)==106
    assert sum(x["source_migration"]==MOD.GROUPS[1] for x in rows)==87
    assert sum(x["n"] for x in rows)==3136
    assert sum(x["far500"] for x in rows)==164
    assert all(0<=x["far500"]<=x["far100"]<=x["n"] for x in rows)
    assert MOD.label(MOD.ROSTER_UNKNOWN)=="NOT_LISTED_UNKNOWN"
    assert MOD.HISTORICAL_CAUSAL_CLAIM is False


def test_published_summary_and_deterministic_bootstraps():
    observed,rows=MOD.analyze()
    expected=json.loads((D/"migrant_roster_host_gap_results.json").read_text(encoding="utf8"))
    assert len(rows)==expected["gap_comparison_species"]==193
    assert observed["full_candidate_butterfly_species"]==expected["full_candidate_butterfly_species"]==225
    assert observed["status"]==expected["status"]
    assert observed["not_listed_is_unknown"] is True
    assert observed["plant_records_can_postdate_butterfly"] is True
    for g in MOD.GROUPS:
        a,b=observed["groups"][g],expected["groups"][g]
        for k in ("species","island_first_index_events","previous_record_gap_500km_events","species_with_500km_gap"):
            assert a[k]==b[k]
        for k in ("gap_500_event_weighted","gap_500_species_equal"):
            close(a[k],b[k])
        for k in a["correlations"]:
            close(a["correlations"][k],b["correlations"][k])
            for actual,published in zip(a["bootstrap_95"][k]["ci95"],b["bootstrap_95"][k]["ci95"]):
                close(actual,published)
            assert a["bootstrap_95"][k]["valid_resamples"]==b["bootstrap_95"][k]["valid_resamples"]==1200
        for cutoff in ("5","10","20"):
            aa,bb=a["minimum_record_thresholds"][cutoff],b["minimum_record_thresholds"][cutoff]
            assert aa["species"]==bb["species"]
            for k in aa["correlations"]:
                close(aa["correlations"][k],bb["correlations"][k])
        aa=observed["followup_retrospective_2001_2024_plant_pooled"][g]
        bb=expected["followup_retrospective_2001_2024_plant_pooled"][g]
        for name in ("ALIEN_CANDIDATE_RECORDED","NATIVE_ONLY_CLASSIFIED"):
            assert aa["host_groups"][name]["n"]==bb["host_groups"][name]["n"]
            assert aa["host_groups"][name]["reindexed"]==bb["host_groups"][name]["reindexed"]
            close(aa["host_groups"][name]["rate"],bb["host_groups"][name]["rate"])
        aaa,bbb=aa["within_species"],bb["within_species"]
        assert aaa["paired_species"]==bbb["paired_species"]
        for k in ("mean_alien_minus_native_reindex_rate",):
            close(aaa[k],bbb[k])
        for x,y in zip(aaa["ci95_resampled"],bbb["ci95_resampled"]):
            close(x,y)


def test_source_joined_193_row_csv():
    expected,_=MOD.join()
    with (D/"migrant_roster_host_gap_193_species.csv").open(newline="",encoding="utf8") as f:
        observed=list(csv.DictReader(f))
    assert len(observed)==len(expected)==193
    for a,b in zip(observed,expected):
        assert a["speciesKey"]==b["speciesKey"]
        assert a["source_migration"]==b["source_migration"]
        assert int(a["n"])==b["n"] and int(a["far500"])==b["far500"]
        close(float(a["host"]),b["host"])
        close(float(a["ports"]),b["ports"])

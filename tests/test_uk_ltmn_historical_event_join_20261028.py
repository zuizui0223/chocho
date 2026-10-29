"""Source-independent boundary protections: missing source rows are never absent."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_ltmn_historical_event_join_20261028 import key_site,integral,time_key,parse_event,joined_source

def visit(site=5,year=1988,day=8,month=6,clock=.5):
    return dict(UKBMS_SITENO=site,DAY=day,MONTH=month,YEAR=year,START_TIME=clock)

def count(site=5,year=1988,day=8,month=6,clock=.5,sp="Meadow brown",n=2):
    return dict(**visit(site,year,day,month,clock),WEEKNO=9,COMMON_NAME=sp,VISIT_TOTAL=n)

def test_source_event_keys_robust_and_date():
    assert parse_event(visit())==("5",1988,6,8,"720")
    assert parse_event(visit(clock="12:00"))==parse_event(visit())
    assert parse_event(visit(month=2,day=30)) is None

def test_actual_zero_must_be_explicit_not_unrecorded_species():
    out=joined_source([visit(),visit(year=2018)],
        [count(sp="Meadow brown",n=3),count(sp="Small white",n=0),count(year=2018,sp="Meadow brown",n=7)])
    q=out["source_counts_qc"]
    assert q["explicit_species_count_zero_records"]==1
    assert q["unique_named_species_count_records"]==3
    assert q["exact_matched_species_count_records_to_actual_visits"]==3
    assert out["inference_firewall"]["missing_named_species_rows_are_not_nondetections"]

def test_mismatched_source_day_is_not_assumed_visit():
    out=joined_source([visit(day=9)],[count(day=8,sp="Meadow brown",n=0)])
    assert out["source_counts_qc"]["source_count_rows_without_visit_on_same_calendar_day"]==1
    assert out["source_counts_qc"]["exact_matched_species_count_records_to_actual_visits"]==0

def test_same_visit_duplicate_count_not_double_counted():
    out=joined_source([visit()],[count(n=4),count(n=4),count(n=6)])
    assert out["source_counts_qc"]["unique_named_species_count_records"]==1
    assert out["source_counts_qc"]["duplicate_identical_count_rows"]==1
    assert out["source_counts_qc"]["nonidentical_duplicate_count_source_rows"]==1

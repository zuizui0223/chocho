"""The deficient-cell audit must remain baseline-only and preserve missingness."""
import importlib.util
from pathlib import Path
import sys

PATH = Path(__file__).parents[1]/"scripts/audit_unseen_island_baseline_20261029.py"
sys.path.insert(0, str(PATH.parent))


def module():
    assert PATH.exists(), "Unseen baseline audit is not implemented"
    spec = importlib.util.spec_from_file_location("baseline_audit", PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def candidate(key="new", **changes):
    return {"island_key": key, "area_km2": 100, "representative_lat": 0,
            "representative_lon": 0, "antarctica_could_be_nearer": 0,
            "mainland_coast_gap_upper_km": 99.9,
            "nearest_wpi_coast_distance_lower_km": 100.1, **changes}


def test_roster_excludes_known_late_outcomes_and_uncertain_boundaries():
    m = module()
    rows = [candidate("new"), candidate("known"), candidate("large", area_km2=20000),
            candidate("small", area_km2=99.9), candidate("gap", mainland_coast_gap_upper_km=100),
            candidate("port", nearest_wpi_coast_distance_lower_km=100),
            candidate("antarctic", antarctica_could_be_nearer=1), candidate("polar", representative_lat=55.01)]
    assert [r["island_key"] for r in m.select_candidates(rows, {"known"})] == ["new"]
    assert m.select_candidates(list(reversed(rows)), {"known"}) == m.select_candidates(rows, {"known"})


def source_result(count=20, assigned=20, years=3, status="COMPLETE"):
    return {"status": status, "count": count,
            "facets": {"SPECIES_KEY": {"A": assigned} if assigned else {},
                       "YEAR": {str(2010+i): 1 for i in range(years)},
                       "DATASET_KEY": {"D": count} if count else {}},
            "unassigned_facet_records": {"SPECIES_KEY": count-assigned}}


def test_baseline_quality_keeps_true_zero_separate_from_failed_query():
    m = module()
    assert m.baseline_quality(source_result())["eligible"]
    assert not m.baseline_quality(source_result(20,17))["eligible"]
    assert not m.baseline_quality(source_result(20,20,2))["eligible"]
    zero = m.baseline_quality(source_result(0,0,0))
    failed = m.baseline_quality(source_result(status="MISSING_QUERY_FAILED"))
    assert zero["records"] == 0
    assert failed["records"] is None
    assert failed["eligible"] is None


def test_audit_cannot_declare_futility_from_incomplete_queries():
    m = module()
    complete_zero = m.baseline_quality(source_result(0,0,0))
    missing = m.baseline_quality(source_result(status="MISSING_QUERY_FAILED"))
    assert m.decision([complete_zero, missing], 2)["status"] == "INCOMPLETE_BASELINE_AUDIT"
    assert m.decision([complete_zero]*2, 2)["status"] == "NO_LATE_ACQUISITION_BASELINE_CELL_INSUFFICIENT"
    good = m.baseline_quality(source_result())
    assert m.decision([good]*3, 3)["status"] == "BASELINE_CELL_SCREEN_CLEARED_COHORT_NOT_YET_FIXED"


def test_query_cannot_silently_acquire_late_outcomes():
    m = module()
    q = m.baseline_payload("POLYGON((0 0,1 0,1 1,0 0))")
    predicates = q["predicate"]["predicates"]
    years = [r for r in predicates if r.get("key") == "YEAR"]
    assert years == [{"type":"greaterThanOrEquals","key":"YEAR","value":"2010"},
                     {"type":"lessThanOrEquals","key":"YEAR","value":"2017"}]
    assert q["limit"] == 0
    assert q["parameters"]["CHECKLIST_KEY"] == ["7ddf754f-d193-4cc9-b351-99906754a03b"]

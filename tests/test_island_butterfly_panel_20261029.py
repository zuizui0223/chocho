"""Failures here protect geographic selection and complete GBIF counts."""
import importlib.util
import sys
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1]/"scripts/acquire_island_butterfly_panel_20261029.py"


def module():
    assert PATH.exists(), "The bounded island acquisition implementation is missing"
    spec = importlib.util.spec_from_file_location("island_butterfly_panel", PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def candidate(key, area, lon):
    return dict(island_key=key, area_km2=area, representative_lon=lon,
                representative_lat=10, antarctica_could_be_nearer=0,
                mainland_coast_gap_lower_km=50, mainland_coast_gap_upper_km=51,
                nearest_wpi_coast_distance_lower_km=10,
                nearest_wpi_coast_distance_upper_km=11)


def response(count, names):
    return dict(offset=0, limit=0, endOfRecords=False, count=count, results=[],
                facets=[dict(field="SPECIES_KEY", counts=[dict(name=n, count=c) for n,c in names])])


def test_geography_selection_is_order_independent_and_second_is_farthest():
    m = module()
    rows = [candidate("a", 90, 0), candidate("b", 80, 1), candidate("c", 70, 90)]
    a = m.select_panel(rows)
    b = m.select_panel(list(reversed(rows)))
    assert [r["island_key"] for r in a] == ["a", "c"]
    assert a == b
    rows[0]["mainland_coast_gap_upper_km"] = 101
    assert "a" not in [r["island_key"] for r in m.select_panel(rows)]


def test_facets_are_paged_to_exhaustion_not_occurrence_endofrecords():
    m = module()
    pages = {0: response(6, [("a",3),("b",2)]), 2: response(6, [("c",1)])}
    answer = m.collect_facets(lambda offset, recheck=False: pages[offset], fields=("SPECIES_KEY",), page_size=2)
    assert answer["count"] == 6
    assert answer["facets"]["SPECIES_KEY"] == {"a":3,"b":2,"c":1}


@pytest.mark.parametrize("pages", [
    {0:response(6,[("a",3),("b",2)]),2:response(7,[("c",1)])},
    {0:response(6,[("a",3),("b",2)]),2:response(6,[("a",1)])},
])
def test_changed_total_and_duplicated_facet_ids_fail_closed(pages):
    m = module()
    with pytest.raises(ValueError):
        m.collect_facets(lambda offset,recheck=False: pages[offset],fields=("SPECIES_KEY",),page_size=2)


def test_zero_requires_successful_response_and_failure_remains_missing():
    m = module()
    zero = m.query_safely(lambda: m.collect_facets(lambda offset,recheck=False:response(0,[]),fields=("SPECIES_KEY",)))
    assert zero["status"] == "COMPLETE" and zero["count"] == 0
    def failed():
        raise TimeoutError("read incomplete")
    bad = m.query_safely(failed)
    assert bad["status"] == "MISSING_QUERY_FAILED" and bad["count"] is None


def test_first_page_movement_is_rejected():
    m = module()
    with pytest.raises(ValueError, match="changed"):
        m.collect_facets(lambda offset,recheck=False:response(2,[("a",1),("b",1)]) if recheck else response(2,[("a",2)]),fields=("SPECIES_KEY",))


def test_lost_dataset_margin_fails_instead_of_certifying_incomplete_facets():
    m = module()
    page = dict(count=5, results=[], facets=[dict(field="DATASET_KEY",counts=[dict(name="d",count=4)])])
    with pytest.raises(ValueError, match="margin"):
        m.collect_facets(lambda offset,recheck=False:page, fields=("DATASET_KEY",))


def test_coverage_keeps_higher_rank_only_records_separate_from_zero_and_failure():
    m = module()
    def row(key,window,count,species):
        return dict(island_key=key,window=window,count=count,status="COMPLETE" if count is not None else "MISSING_QUERY_FAILED",
                    facets={"SPECIES_KEY":species},unassigned_facet_records={"SPECIES_KEY":count-sum(species.values()) if count is not None else 0})
    data = [row("a","2010_2017",5,{"1":2,"2":1,"3":1}),row("a","2018_2025",3,{"1":1,"2":1,"3":1}),
            row("b","2010_2017",2,{}),row("b","2018_2025",0,{}),row("c","2010_2017",None,{}),row("c","2018_2025",0,{})]
    summary = m.coverage_summary(data,[{"island_key":k} for k in ("a","b","c")])
    assert summary["islands_atleast3_species_bins_each_window"] == 1
    assert summary["positive_records_without_any_species_bin_windows"] == 1
    assert summary["complete_zero_record_island_windows"] == 2
    assert summary["failed_island_windows"] == 1
    assert summary["records_without_species_bin"] == 3
